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

argument_parser = argparse.ArgumentParser(
        description='Run foam experiment of increasing size and record the time took to execute the'
            'Fast Multipole Method')
argument_parser.add_argument('lattice_max_size', type=int)
args = argument_parser.parse_args()

sizes = [ i for i in range(1, args.lattice_max_size + 1) ]

config_file = pathlib.Path('assets/bubblelattice.txt')
config = config_file.read_text()

working_directory = pathlib.Path('BubbleLattice')
working_directory.mkdir(exist_ok=True)

for size in sizes:
    experiment_path = working_directory / f'Size{size}'
    shutil.rmtree(experiment_path.as_posix())
    experiment_path.mkdir()

    experiment_config_file = experiment_path / 'config.txt'
    experiment_config_file.write_text(config.format(bubble_lattice_size=size))

    completed_process = subprocess.run(
            ['./SoapFilm3D', experiment_config_file.as_posix(), 'headless'],
            capture_output=True,
            text=True)

    getMostRecentOutputDirectory().rename(experiment_path / 'output')

    if completed_process.returncode != 0:
        (experiment_path / 'error').touch()
        continue

    fmm_execution_time = []
    for line in completed_process.stdout.split('\n'):
        if line.startswith('FMMExecution'):
            fmm_execution_time.append(float(line.split()[-1]))

    (experiment_path / 'result.csv').write_text(
        '\n'.join(map(str, fmm_execution_time)))

