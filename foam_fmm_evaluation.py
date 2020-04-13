import pathlib
import subprocess
import argparse
import shutil
import csv

CONFIG_FILE_PATH = pathlib.Path('assets/bubblelattice.txt')
RESULT_FILE_NAME = 'result.csv'

def getOutputDirectoryTime(output_directory_path):
    return int(output_directory_path.name[len('output_'):])

def getMostRecentOutputDirectory():
    most_recent = None
    for cwd_content in pathlib.Path.cwd().iterdir():
        if not cwd_content.is_dir():
            continue

        if not cwd_content.name.startswith('output_'):
            continue

        if (
                most_recent is None
                or getOutputDirectoryTime(most_recent)
                    < getOutputDirectoryTime(cwd_content)):
            most_recent = cwd_content

    return most_recent

def getTimePrefix(method):
    if method == 'fmmtl':
        return 'FMMExecution'
    else:
        return 'NaiveExecution'

def getMethodSubdirectory(method):
    if method == 'fmmtl':
        return 'FMM'
    else:
        return 'Naive'

def getDataPrefixes(method):
    return [ getTimePrefix(method), 'NumberVertices', 'VertexDensity' ]

argument_parser = argparse.ArgumentParser(
        description='Run foam experiment of increasing size and record the time took to execute the'
            'Fast Multipole Method')
argument_parser.add_argument('lattice_max_size', type=int)
argument_parser.add_argument('-m', '--method', action='append', choices=['fmmtl', 'naive'])
argument_parser.add_argument('-s', '--save-mesh', action='store_true')
argument_parser.add_argument('-T', '--save-mesh-period', type=int, default=1)
argument_parser.add_argument('-o', '--output-directory', default='BubbleLattice')
argument_parser.add_argument('-c', '--config', default='assets/bubblelattice.txt')
args = argument_parser.parse_args()

sizes = [ i for i in range(1, args.lattice_max_size + 1) ]

config = pathlib.Path(args.config).read_text()
output_directory = pathlib.Path(args.output_directory)

output_directory.mkdir(parents=True,exist_ok=True)

for size in sizes:
    experiment_size_path = experiment_path = output_directory / f'Size{size}'
    for method in args.method:
        print(f'Starting experiment of size {size}.')
        experiment_path = experiment_size_path / getMethodSubdirectory(method)
        if experiment_path.exists():
            print(f'Output directory {experiment_path.as_posix()} was present, deleting.')
            shutil.rmtree(experiment_path.as_posix())
        experiment_path.mkdir(parents=True)

        experiment_config_file = experiment_path / 'config.txt'
        print(f'Saving configuration at {experiment_config_file.as_posix()}')
        experiment_config_file.write_text(config.format(
            bubble_lattice_size=size,
            fmmtl_enable=int(method == 'fmmtl'),
            output_mesh=int(args.save_mesh),
            output_mesh_every_n_frames=args.save_mesh_period))

        print(f'Starting the simulation')
        completed_process = subprocess.run(
                ['./SoapFilm3D', experiment_config_file.as_posix(), 'headless'],
                capture_output=True,
                text=True)

        output_path = getMostRecentOutputDirectory()
        if args.save_mesh:
            output_path.rename(experiment_path / 'output')
        else:
            shutil.rmtree(output_path.as_posix(), ignore_errors=True)
            

        if completed_process.returncode != 0:
            print('An error as occured during the simulation, skiping')
            (experiment_path / 'error').touch()
            continue

        print('Writing result file')
        data_prefixes = getDataPrefixes(method)
        data = { data_prefix : [] for data_prefix in data_prefixes }
        for line in completed_process.stdout.split('\n'):
            line_components = line.split()
            if len(line_components) != 2:
                continue
            prefix = line_components[0]
            if prefix in data_prefixes:
                data[prefix].append(float(line_components[-1]))

        with open((experiment_path / 'result.csv').as_posix(), mode='w') as result_file:
            result_object = csv.writer(result_file, delimiter=' ')
            result_object.writerow(getDataPrefixes(method))
            result_object.writerows(zip(*[ data[prefix] for prefix in data_prefixes]))

