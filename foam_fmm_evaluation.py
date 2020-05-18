#!/usr/bin/python3.8

import pathlib
import subprocess
import argparse
import shutil
import csv
import itertools
import functools

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



class SoapFilmSimulation(object):

    config_file_name = 'config.txt'

    def __init__(self,
            output_directory,
            delete_existing=False,
            config=None,
            headless=True,
            capture_stdout=True,
            timeout=None
            ):
        self.capture_stdout = capture_stdout
        self.headless = headless
        self.timeout = timeout

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

    @property
    def config_file(self):
        return self.output_directory / SoapFilmSimulation.config_file_name

    def writeConfig(self):
        self.config.writeToFile(self.config_file.as_posix())

    def run(self):
        completed_process = subprocess.run(
                [
                    './SoapFilm3D',
                    self.config_file.as_posix(),
                    'headless' if self.headless else 'output'
                ],
                capture_output=self.capture_stdout,
                text=True,
                timeout=self.timeout)

        if completed_process.returncode != 0:
            print(completed_process.stderr)
            raise RuntimeError('An error as occured during the simulation')

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
        return [ 
                key
                for key in self.options
                if len(self.options[key]) > 1 and key not in self.order
            ]

    def __iter__(self):
        return SimulationParameterProductIterator(self.options, self.getRelevantConfigurationKeys())

    def __repr__(self):
        return repr(self.options)

    def __str__(self):
        return str(self.options)

if  __name__ == '__main__':

    argument_parser = argparse.ArgumentParser(
            description='Run foam experiment of increasing size and record the time took to execute'
                'the Fast Multipole Method')
    argument_parser.add_argument(
            '-m', '--method', action='append', choices=['fmmtl', 'naive'], default=[])
    argument_parser.add_argument('-s', '--save-mesh', action='store_true')
    argument_parser.add_argument('-T', '--save-mesh-period', type=int, default=1)
    argument_parser.add_argument('-c', '--config')
    argument_parser.add_argument('-t', '--simulation-time', type=float, default=4.0)
    argument_parser.add_argument('-o', '--sim-option', action='append', default=[])
    argument_parser.add_argument('--no-save-stdout' , action='store_true')
    argument_parser.add_argument('--no-run', action='store_true',
        help='The simulaton is not run, only the configuration file is saved.')
    argument_parser.add_argument('--no-headless', action='store_true')
    argument_parser.add_argument('--directory-order', action='append', default=[])
    argument_parser.add_argument('--load',
        help='Load simulations produced with this program. Implies -o scene=load --no-headless'
            ' --no-save-stdout. The options must be specified in the same order for this to work' 
            ' correctly.')
    argument_parser.add_argument('--timeout', type=float,
            help='Stops each simulation after TIMEOUT seconds.')
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

    config.output_mesh = args.save_mesh
    config.output_mesh_every_n_frames = args.save_mesh_period
    config.load_increment = args.save_mesh_period
    config.simulation_time = args.simulation_time

    simulation_parameter_product = SimulationParameterProduct(args.directory_order)
    for sim_option in args.sim_option:
        key, values = sim_option.split('=')
        values = values.split(',')
        simulation_parameter_product.addOptions(key, values)

    for path, options in simulation_parameter_product:
        print(f'Starting simulation with options {options}.')
        experiment_path = pathlib.Path(args.output_directory) / path

        if args.load is not None:
            config.load_dir = pathlib.Path(args.load) / path / 'output'
        config.update(options)

        simulation = SoapFilmSimulation(
                experiment_path,
                delete_existing=True,
                capture_stdout=not args.no_save_stdout,
                config=config,
                headless=not args.no_headless,
                timeout=args.timeout)

        if args.no_run:
            continue

        try:
            stdout = simulation.run()
        except RuntimeError as e:
            (experiment_path / 'error').touch()
            print(e)
            continue
        except subprocess.TimeoutExpired as timeout_expired:
            stdout = timeout_expired.stdout.decode()

        if not args.no_save_stdout:
            (experiment_path / 'stdout').write_text(stdout)

