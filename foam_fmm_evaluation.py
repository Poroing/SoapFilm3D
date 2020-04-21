import pathlib
import subprocess
import argparse
import shutil
import csv
import itertools

class SoapFilmSimulationConfigFile(object):

    TRANSLATION = {
            'subdivisions' : 'mesh-size-n',
            'size' : 'mesh-size-m'
            }

    @staticmethod
    def getTranslation(key):
        return SoapFilmSimulationConfigFile.TRANSLATION.get(
                key,
                key)

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

    def __setattr__(self, name, value):
        if name == 'config':
            object.__setattr__(self, name, value)
        else:
            self.__setitem__(self.getConfigKeyFromVariableName(name), value)

    def __getattr__(self, name):
        return self.__getitem__(self.getConfigKeyFromVariableName(name))

    def __getitem__(self, key):
        return self.config[SoapFilmSimulationConfigFile.getTranslation(key)]

    def __setitem__(self, key, value):
        self.config[SoapFilmSimulationConfigFile.getTranslation(key)] = value

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

    def __init__(self, output_directory, delete_existing=False, config=None):
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
                ['./SoapFilm3D', self.config_file.as_posix(), 'headless'],
                capture_output=True,
                text=True)

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

if  __name__ == '__main__':


    argument_parser = argparse.ArgumentParser(
            description='Run foam experiment of increasing size and record the time took to execute'
                'the Fast Multipole Method')
    argument_parser.add_argument(
            '-m', '--method', action='append', choices=['fmmtl', 'naive'], default=[])
    argument_parser.add_argument('-s', '--save-mesh', action='store_true')
    argument_parser.add_argument('-T', '--save-mesh-period', type=int, default=1)
    argument_parser.add_argument('-o', '--output-directory', default='BubbleLattice')
    argument_parser.add_argument('-c', '--config')
    argument_parser.add_argument('-r', '--resolution', action='append', type=float, default=[])
    argument_parser.add_argument('-S', '--size', action='append', type=int, default=[])
    argument_parser.add_argument('-t', '--simulation-time', type=float, default=4.0)
    argument_parser.add_argument('--scene', default='bubblelattice')
    argument_parser.add_argument('--subdivisions', type=int, default=2)
    args = argument_parser.parse_args()

    if args.method == []:
        args.method = ['fmmtl']

    if args.config is not None:
        config = SoapFilmSimulationConfigFile.fromConfigFile(args.config)
    else:
        config = SoapFilmSimulationConfigFile()

    config.scene = args.scene
    config.mesh_size_n = args.subdivisions
    config.output_mesh = args.save_mesh
    config.output_mesh_every_n_frames=args.save_mesh_period
    config.simulation_time = args.simulation_time

    for size, resolution, method in itertools.product(args.size, args.resolution, args.method):
        print(f'Starting experiment of size {size}, resolution {resolution} and method {method}.')
        experiment_path = (
                output_directory
                / f'Resolution{resolution}'
                / f'Size{size}'
                / method.capitalize()
                )

        config.mesh_size_m = size
        config.remeshing_resolution = resolution
        config.fmmtl = (method == 'fmmtl')

        simulation = SoapFilmSimulation(experiment_path, delete_existing=True, config=config)
        stdout = simulation.run()
        (experiment_path / 'stdout').write_text(stdout)

