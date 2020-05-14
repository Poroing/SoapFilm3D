from foam_fmm_evaluation import *
import pathlib
import argparse
from divide_frame_number import dividesFramesNumber 
import ffmpeg
import csv

def checkPathExists(path):
    if not path.exists():
        raise ValueError(f'{path.as_posix()} does not exists.')

def video(path, config, skip_existing=False, overwrite=False):
    video_path = path / 'video.webm'
    if video_path.exists() and skip_existing:
        return

    temporary_frame_directory = path / 'Frames'
    divider = int(config.get('output-png-every-n-frames',
                config.get('output-mesh-every-n-frames',
                    1)))

    dividesFramesNumber(
            'frame',
            '.png',
            temporary_frame_directory,
            path / 'output',
            divider)

    frame_rate = int(1 / (divider * config.get('time-step', 0.01)))

    stream = ffmpeg.input((temporary_frame_directory / 'frame%d.png').as_posix(), r=frame_rate)
    stream = ffmpeg.output(stream, video_path.as_posix())
    if overwrite:
        stream = ffmpeg.overwrite_output(stream)
    ffmpeg.run(stream)

def getObjs(path):
    output_path = path / 'output'
    for file in output_path.iterdir():
        if file.suffix == '.obj':
            yield file

def getFrameFromObj(path):
    return path.parent / f'frame{int(path.name[4:10])}.png'


def render(path, skip_existing=False, show_blender_stdout=False, number_subdivisions=0,scale=1):
    for file in getObjs(path):

        output_frame = getFrameFromObj(file)
        if skip_existing and output_frame.exists():
            continue

        print(f'{file.as_posix()} --> {output_frame.as_posix()}')
        subprocess.run(
                [ 
                    'primusrun',
                    'blender',
                    '--background',
                    '--python',
                    'render.py',
                    '--',
                    output_frame.as_posix(),
                    file.as_posix(),
                    '--number-subdivisions', str(number_subdivisions),
                    '--scale', str(scale)
                ],
                capture_output=not show_blender_stdout
            )

def deleteRender(path):
    for file in getObjs(path):
        getFrameFromObj(file).unlink(missing_ok=True)

def produceCSV(path, skip_existing=False, delimiter=' '):
    csv_path = path / 'data.csv'
    if csv_path.exists() and skip_existing:
        return

    TOKEN = ['NumberVertices', 'FMMExecution', 'NaiveExecution']
    data = {}
    stdout = (path / 'stdout').read_text()
    for line in stdout.split('\n'):
        line_tokens = line.split(' ')
        if line_tokens[0] in TOKEN:
            data.setdefault(line_tokens[0], []).append(float(line_tokens[1]))

    data_length = None
    for value in data.values():
        if data_length is None:
            data_length = len(value)

        if data_length != len(value):
            raise ValueError('The standard output does not contain as much data of each type')

    data_keys = data.keys()
    data = [ { k : v[i] for k, v in data.items() } for i in range(data_length) ]
    with open(csv_path.as_posix(), mode='w', newline='') as csv_file:
        csv_writer = csv.DictWriter(csv_file, fieldnames=data_keys, delimiter=args.csv_delimiter)
        csv_writer.writeheader()
        csv_writer.writerows(data)

def getDescriptionText(config, relevant_configuration_keys):
    return '\n'.join([ f'{key} : {config[key]}' for key in relevant_configuration_keys ])


def concatenateVideos(
        paths,
        configs,
        relevant_configuration_keys,
        output_path,
        overwrite=False,
        font_color='black'):

    streams = []
    for path, config in zip(paths, configs):
        stream = ffmpeg.input((path / 'video.webm').as_posix())
        stream = ffmpeg.drawtext(
                stream,
                text=getDescriptionText(config, relevant_configuration_keys),
                fontcolor=font_color)
        streams.append(stream)

    
    stream = ffmpeg.concat(*streams)
    if overwrite:
        stream = ffmpeg.overwrite_output(stream)
    stream = ffmpeg.output(stream, output_path.as_posix())
    ffmpeg.run(stream)

argument_parser = argparse.ArgumentParser()
argument_parser.add_argument('-o', '--sim-option', action='append', default=[])
argument_parser.add_argument('--skip-existing', action='store_true',
    help='Takes precedence over --overwrite')
argument_parser.add_argument('--overwrite', action='store_true')
argument_parser.add_argument('--directory-order', action='append', default=[])
argument_parser.add_argument('--show-blender-stdout', action='store_true')
argument_parser.add_argument('--font-color', default='black')
argument_parser.add_argument('command', choices=['video', 'render', 'delete-render', 'csv'])
argument_parser.add_argument('--number-subdivisions-render', type=int, default=0)
argument_parser.add_argument('--csv-delimiter', default=' ')
argument_parser.add_argument('--scale', type=float, default=1)
argument_parser.add_argument('output_directory')
args = argument_parser.parse_args()

output_directory = pathlib.Path(args.output_directory)
checkPathExists(output_directory)

simulation_parameter_product = SimulationParameterProduct(args.directory_order)
for sim_option in args.sim_option:
    key, values = sim_option.split('=')
    values = values.split(',')
    simulation_parameter_product.addOptions(key, values)

paths = []
configs = []
for path, _ in simulation_parameter_product:
    experiment_path = output_directory / path
    checkPathExists(experiment_path)
    paths.append(experiment_path)

    config_path = experiment_path / 'config.txt'
    checkPathExists(config_path)
    config = SoapFilmSimulationConfigFile.fromConfigFile(config_path.as_posix())
    configs.append(config)

    if args.command == 'video':
        video(experiment_path, config, skip_existing=args.skip_existing, overwrite=args.overwrite)

    if args.command == 'render':
        render(
                experiment_path,
                skip_existing=args.skip_existing,
                show_blender_stdout=args.show_blender_stdout,
                number_subdivisions=args.number_subdivisions_render,
                scale=args.scale)

    if args.command == 'delete-render':
        deleteRender(experiment_path)

    if args.command == 'csv':
        produceCSV(experiment_path, skip_existing=args.skip_existing, delimiter=args.csv_delimiter)

if args.command == 'video' and simulation_parameter_product.getNumberDifferentConfigurations() > 1:
    concatenateVideos(
            paths,
            configs,
            simulation_parameter_product.getRelevantConfigurationKeys(),
            output_directory / 'video.webm',
            font_color=args.font_color)

        
