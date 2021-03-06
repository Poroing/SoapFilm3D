#!/usr/bin/python3

import pathlib
import subprocess
import argparse
import shutil
import csv
import itertools
import functools
import time
import collections
import os
import collections.abc
import textwrap

class FailedSimulation(Exception):

    def __init__(self, stdout, stderr):
        self.stdout = stdout
        self.stderr = stderr

class SoapFilmSimulationConfigFile(object):

    @staticmethod
    def getConfigKeyFromVariableName(variable_name):
        return variable_name.replace('_', '-')

    @staticmethod
    def getVariableNameFromConfigKey(config_key):
        return config_key.replace('-', '_')

    @classmethod
    def fromConfigFile(cls, config_file):
        with open(config_file, newline='') as config_file:
            config_object = csv.reader(config_file, delimiter=' ')
            a = dict(config_object)
            return cls(**a)

    def __init__(self, **kwargs):
        self.config = { self.getConfigKeyFromVariableName(k) : v for k, v in kwargs.items() }

    def update(self, options):
        self.config.update(options)

    def get(self, key, default):
        return self.config.get(key, default)

    def __setattr__(self, name, value):
        if name == 'config':
            object.__setattr__(self, name, value)
        else:
            self.__setitem__(self.getConfigKeyFromVariableName(name), value)

    def __getattr__(self, name):
        return self.__getitem__(self.getConfigKeyFromVariableName(name))

    def __getitem__(self, key):
        return self.config[key]

    def __setitem__(self, key, value):
        self.config[key] = value

    def __getstate__(self):
        """ Make sure that pickle module use __dict__ for picking."""
        return self.__dict__

    def __setstate__(self, state):
        self.__dict__.update(state)


    def writeToFile(self, filename):
        with open(filename, mode='w', newline='') as config_file:
            config_object = csv.writer(config_file, delimiter=' ')
            for key, value in self.config.items():
                if isinstance(value, bool):
                    config_object.writerow((key, str(int(value))))
                else:
                    config_object.writerow((key, str(value)))

def isGlobalTimeoutReached(end_time):
    return end_time is not None and time.time() > end_time



class SoapFilmSimulation(object):

    config_file_name = 'config.txt'

    @staticmethod
    def hasFinishedInit(stdout):
        return stdout.find('Execution') != -1

    @staticmethod
    def hasNotEnoughSphere(stderr):
        return stderr.find('Not enough spheres') != -1

    def __init__(self,
            output_directory,
            simulation_executable,
            delete_existing=False,
            config=None,
            headless=True,
            capture_stdout=True,
            timeout=None,
            retry_on_failed_initialization=False,
            retry_on_not_enough_spheres=False,
            profile=False
            ):
        self.simulation_executable = simulation_executable
        self.capture_stdout = capture_stdout
        self.headless = headless
        self.timeout = timeout
        self.retry_on_failed_initialization = retry_on_failed_initialization
        self.retry_on_not_enough_spheres = retry_on_not_enough_spheres
        self.profile = profile

        if config is None:
            config = SoapFilmSimulationConfigFile()
        self.config = config

        self.output_directory = pathlib.Path(output_directory)
        if self.output_directory.exists():
            if delete_existing:
                shutil.rmtree(self.output_directory.as_posix())
            else:
                raise ValueError(f'The output directory {self.output_directory.as_posix()} exists')
        self.output_directory.mkdir(parents=True)

        self.config.output_dir = output_directory / 'output'
        self.writeConfig()

    def shouldCaptureStdout(self):
        return (
                self.capture_stdout
                or self.retry_on_failed_initialization
                or self.retry_on_not_enough_spheres
            )


    @property
    def config_file(self):
        return self.output_directory / SoapFilmSimulation.config_file_name

    def writeConfig(self):
        self.config.writeToFile(self.config_file.as_posix())

    def shouldRetry(self, stdout, stderr):
        return (
            (
                self.retry_on_failed_initialization
                and not SoapFilmSimulation.hasFinishedInit(stdout)
            )
            or
            (
                self.retry_on_not_enough_spheres
                and SoapFilmSimulation.hasNotEnoughSphere(stderr)
            )
        )


    def run(self):
        number_tries = 1
        environment_variable = {}
        environment_variable.update(os.environ)
        if self.profile:
            environment_variable['CPUPROFILE'] = str(self.output_directory / 'prof.out')
            environment_variable['OMP_NUM_THREADS'] = str(1)
        while True:
            print(f'Try {number_tries}')
            completed_process = subprocess.run(
                    [
                        self.simulation_executable,
                        self.config_file.as_posix(),
                        'headless' if self.headless else 'output'
                    ],
                    capture_output=self.shouldCaptureStdout(),
                    text=True,
                    timeout=self.timeout,
                    env=environment_variable)

            if completed_process.returncode != 0:
                if not SoapFilmSimulation.hasFinishedInit(completed_process.stdout):
                    print('Did Not Finished Initialization')
                if SoapFilmSimulation.hasNotEnoughSphere(completed_process.stderr):
                    print('Not enough spheres')

                if self.shouldRetry(completed_process.stdout, completed_process.stderr):
                    number_tries += 1
                    continue
                else:
                    print(completed_process.stderr)
                    raise FailedSimulation(completed_process.stdout, completed_process.stderr)
            break

        return completed_process.stdout

class SoapFilmSimulationResult(object):

    def __init__(self, output_directory, stdout=None):
        self.output_directory = pathlib.Path(output_directory)
        if not self.output_directory.exists():
            raise ValueError(
                    f'Output directory f{self.output_directory.as_posix()} does not exists')
        
        config_file = self.output_directory / SoapFilmSimulation.config_file_name
        if not config_file.exists:
            raise ValueError(f'Config file f{config_file.as_posix()} does not exists')
        self.config = SoapFilmSimulationConfigFile.fromConfigFile(config_file.as_posix())

        if stdout is None:
            stdout = (self.output_directory / 'stdout').read_text()
        self.stdout = stdout

class SimulationParameterProductIterator(object):

    def __init__(self, options, order):
        self.order = order
        self.options = options
        self.iterator = itertools.product(
                    *[ [ (key, value) for value in self.options[key] ] for key in self.options if
                        len(self.options[key]) > 0 ]
                )

    def __iter__(self):
        return self

    def __next__(self):
        next_options = dict(next(self.iterator))
        path = pathlib.Path('')
        for key in self.order:
            path = path / (key.capitalize() + str(next_options[key]).capitalize())
        return path, next_options


class SimulationParameterProduct(object):

    def __init__(self, order=None):
        if order is None:
            order = []
        self.order = order
        self.options = {}

    def getNumberDifferentConfigurations(self):
        return functools.reduce(
                lambda x, y: x * y,
                map(len, self.options.values()),
                1)

    def addOption(self, key, value):
        self.options.setdefault(key, []).append(value)

    def addOptions(self, key, values):
        self.options.setdefault(key, []).extend(values)

    def getRelevantConfigurationKeys(self):
        return self.order + self.getRelevantConfigurationKeysNotInOrder()

    def getRelevantConfigurationKeysNotInOrder(self):
        return self.getRelevantConfigurationKeysNotInList(self.order)

    def getRelevantConfigurationKeysNotInList(self, l):
        return [ 
                key
                for key in self.options
                if len(self.options[key]) > 1 and key not in l
            ]


    def iterateOnSets(self, parameters):
        """Returns a generator on sets of paths. The iteration is on the product of the parameters
            not in options. Each set of paths is represented by a list of pairs, the first element
            of each pair is the path, and the second element is a dictionary that stores the value
            of options related to this path
        """

        #Make sure that the given options are iterated on last.
        options = { key : values for key, values in self.options.items() if key not in parameters }
        options.update(self.options)
        jump_size = functools.reduce(
                lambda v, l: len(l) * v,
                [ self.options[key] for key in parameters ],
                1
            )

        paths_and_configs = list(SimulationParameterProductIterator(
            options,
            self.getRelevantConfigurationKeys()))

        for i in range(0, len(paths_and_configs), jump_size):
            yield paths_and_configs[i : i + jump_size]

    def __iter__(self):
        return SimulationParameterProductIterator(self.options, self.getRelevantConfigurationKeys())

    def __repr__(self):
        return repr(self.options)

    def __str__(self):
        return str(self.options)

if  __name__ == '__main__':

    argument_parser = argparse.ArgumentParser(
            description=textwrap.dedent('''\
                Run soap films simulations on the product of sets of simulation options.
                Example:

                    ./run_simulations.py \\
                        -o mesh-size-n=3 \\
                        -o mesh-size-m=1,2 \\
                        -o scene=mergedbubblelattice,2dbubblelattice

                runs the soap film simulation with options

                    (mesh-size-n=2, mesh-size-m=1, scene=mergedbubblelattice)
                    (mesh-size-n=2, mesh-size-m=1, scene=2dbubblelattice)
                    (mesh-size-n=2, mesh-size-m=2, scene=mergedbubblelattice)
                    (mesh-size-n=2, mesh-size-m=2, scene=2dbubblelattice)

                and stores the standard output, config files, objs, png and rec in respectively:

                    <output-directory>/Mesh-size-m1/SceneMergedbubblelattice
                    <output-directory>/Mesh-size-m1/Scene2dbubblelattice
                    <output-directory>/Mesh-size-m2/SceneMergedBubbleLattice
                    <output-directory>/Mesh-size-m2/Scene2dbubblelattice

                Note that there is no directory Mesh-size-n2 as mesh-size-n only takes one value.
                The standard output is saved in a text file with name stdout.
                The config file is saved in a text file with name config.txt.
                The objs,png and rec are stored as the simulation would do in a subdirectory.
                output.
                '''),
            formatter_class=argparse.RawTextHelpFormatter
        )
    argument_parser.add_argument('-T', '--save-mesh-period', type=int, default=1,
            help='Save mesh every T frames')
    argument_parser.add_argument('-c', '--config',
        help=textwrap.dedent('''\
            Base config file. The other options override the options defined in the given config
            file.
        ''')
        )
    argument_parser.add_argument('-t', '--simulation-time', type=float, default=4.0,
            help=textwrap.dedent('''
                Simulation duration. Note that this is not the timeout but the time used to compute
                how much time-steps are made
            ''')
        )
    argument_parser.add_argument('-o', '--sim-option', action='append', default=[],
            metavar='<keyword>=<value>|<keyword>=<value>,...,<value>|<keyword>=E:<expression>',
            help=textwrap.dedent('''\
                Add given values to the option <keyword>. The first and second form add the given
                values to the option. The third form evaluate the python expression <expression>
                which should evaluate to a list that is added to the option.
            ''')
            )
    argument_parser.add_argument('--no-save-stdout' , action='store_true',
            help=textwrap.dedent('''\
                Don\'t save standard output in a file and print it out directly to standard output
            '''))
    argument_parser.add_argument('--profile', action='store_true',
        help=textwrap.dedent('''\
            Profile the simulation. To use this the executable must have been link with -lprofiler.
            This causes the simulation to stop multithreading and to save the profiling data in the
            related directory as prof.out. See the documentation on gperftools to understand how to
            use this file.
        ''')
        )
    argument_parser.add_argument('--no-run', action='store_true',
        help=textwrap.dedent('''\
            The simulaton is not run, only the configuration file is saved.
        ''')
        )
    argument_parser.add_argument('--no-headless', action='store_true',
        help='A window with the simulation is created')
    argument_parser.add_argument('--directory-order', action='append', default=[],
            help=textwrap.dedent('''\
                Specify in which order the directories should be saved. Each options with more
                than one value must be specified.
                Example:

                    ./run_simulations.py \\
                        -o mesh-size-m=1,2 \\
                        -o scene=mergedbubblelattice,2dbubblelattice\\
                        --directory-order scene --directory-order mesh-size-m

                This will save the results in directories\

                    <output_directory>/Scene<value>/Mesh-size-m<value>

                instead of

                    <output_directory>/Mesh-size-m<value>/Scene<value>/
            ''')
            )
    argument_parser.add_argument('--retry-on-failed-initialization', action='store_true',
        help=textwrap.dedent('''\
            This causes the simulation to be rerun if it failed during the initialization (before it
            has simulated the first step. Implies no --no-save-stdout.
        ''')
        )
    argument_parser.add_argument('--retry-on-not-enough-spheres', action='store_true',
        help=textwrap.dedent('''
            Some scenes, such as newfoam can throw a std::runtime_error because they were not
            able to produce the asked number of sphere. This causes to retry the simulation when
            it happen. Implies no --no-save-stdout
        ''')
        )
    argument_parser.add_argument('--load',
        help=textwrap.dedent('''\
            Load simulations produced with this program. The -T argument should be the same used to
            save the simulation. Implies -o scene=load --no-headless --no-save-stdout. The options
            must be specified in the same order for this to work correctly.
            ''')
        )
    argument_parser.add_argument('--timeout', type=float,
            help='Stops each simulation after TIMEOUT seconds.')
    argument_parser.add_argument('--global-timeout', type=float,
        help='Stop every simulation when timeout GLOBAL_TIMEOUT has passed')
    argument_parser.add_argument('--unstability-algorithm', type=int, metavar='N',
            help=textwrap.dedent('''\
                Starts with the given simulation timeout and double the timeout each time N
                consecutive simulation timeout. This allows to run as many simulation as you want
                without unstable taking all the allocated time.
            ''')
            )
    argument_parser.add_argument(
            '--simulation-executable',
            default=pathlib.Path(__file__).parent / 'build' / 'Apps' / 'SoapFilm3D' / 'SoapFilm3D')
    argument_parser.add_argument('output_directory')
    args = argument_parser.parse_args()

    if args.load is not None:
        args.sim_option.append('scene=load')
        args.no_headless = True
        args.no_save_stdout = True

    if args.config is not None:
        config = SoapFilmSimulationConfigFile.fromConfigFile(args.config)
    else:
        config = SoapFilmSimulationConfigFile()

    if args.unstability_algorithm is not None and args.timeout is None:
        print('--unstability-algorithm needs the --timeout argument')
        exit(1)


    config.output_mesh_every_n_frames = args.save_mesh_period
    config.load_increment = args.save_mesh_period
    config.simulation_time = args.simulation_time

    simulation_parameter_product = SimulationParameterProduct(args.directory_order)
    for sim_option in args.sim_option:
        key, values = sim_option.split('=')
        if values.startswith('E:'):
            values = list(eval(values[2:]))
        else:
            values = values.split(',')
        simulation_parameter_product.addOptions(key, values)

    end_time = None
    if args.global_timeout is not None:
        end_time = time.time() + args.global_timeout

    simulation_timeout = args.timeout

    last_timedout_simulations_parameter = collections.deque()
    simulation_parameter_to_rerun = collections.deque()
    simulation_parameter_iterator = iter(simulation_parameter_product)

    while True:
        if isGlobalTimeoutReached(end_time):
            print('Global timeout attained')
            break

        if len(simulation_parameter_to_rerun) > 0:
            path, options = simulation_parameter_to_rerun.popleft()
        else:
            try:
                path, options = next(simulation_parameter_iterator)
            except StopIteration:
                break

        print(f'Starting simulation with options {options}.')
        experiment_path = pathlib.Path(args.output_directory) / path

        if args.load is not None:
            config.load_dir = pathlib.Path(args.load) / path / 'output'
        config.update(options)

        timeout = None
        if args.timeout is not None and args.global_timeout is not None:
            timeout = min(simulation_timeout, end_time - time.time())
        elif args.timeout is not None:
            timeout = simulation_timeout
        elif args.global_timeout is not None:
            timeout = end_time - time.time()

        simulation = SoapFilmSimulation(
                experiment_path,
                args.simulation_executable,
                delete_existing=True,
                capture_stdout=not args.no_save_stdout,
                config=config,
                headless=not args.no_headless,
                timeout=timeout,
                retry_on_failed_initialization=args.retry_on_failed_initialization,
                retry_on_not_enough_spheres=args.retry_on_not_enough_spheres,
                profile=args.profile)

        if args.no_run:
            continue

        try:
            stdout = simulation.run()
            last_timedout_simulations_parameter.clear()
        except FailedSimulation as e:
            (experiment_path / 'error').touch()
            stdout = e.stdout
            if not args.no_save_stdout:
                stderr = e.stderr
                (experiment_path / 'stderr').write_text(stderr)
        except subprocess.TimeoutExpired as timeout_expired:
            print('Timeout')
            stdout = timeout_expired.stdout.decode()
            if args.unstability_algorithm is not None:
                last_timedout_simulations_parameter.append((path, options))
                if len(last_timedout_simulations_parameter) >= args.unstability_algorithm:
                    simulation_parameter_to_rerun = last_timedout_simulations_parameter.copy()
                    last_timedout_simulations_parameter.clear()
                    simulation_timeout *= 2
                    print(f'Increasing timeout to {simulation_timeout}')

        if not args.no_save_stdout:
            (experiment_path / 'stdout').write_text(stdout)

