import pathlib
import subprocess
import argparse
import shutil

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

argument_parser = argparse.ArgumentParser(
        description='Run foam experiment of increasing size and record the time took to execute the'
            'Fast Multipole Method')
argument_parser.add_argument('lattice_max_size', type=int)
argument_parser.add_argument('-m', '--method', action='append', choices=['fmmtl', 'naive'])
argument_parser.add_argument('-s', '--save-mesh', action='store_true')
args = argument_parser.parse_args()

sizes = [ i for i in range(1, args.lattice_max_size + 1) ]

config_file = pathlib.Path('assets/bubblelattice.txt')
config = config_file.read_text()

working_directory = pathlib.Path('BubbleLattice')
working_directory.mkdir(exist_ok=True)

for size in sizes:
    experiment_size_path = experiment_path = working_directory / f'Size{size}'
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
            output_mesh=int(args.save_mesh)))

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
        execution_time = []
        for line in completed_process.stdout.split('\n'):
            if line.startswith(getTimePrefix(method)):
                execution_time.append(float(line.split()[-1]))

        (experiment_path / 'result.csv').write_text(
            '\n'.join(map(str, execution_time)))

