from foam_fmm_evaluation import *
import pathlib
import argparse
from divide_frame_number import dividesFramesNumber 
import csv
import sys
import multiprocessing

def checkPathExists(path):
    if not path.exists():
        raise ValueError(f'{path.as_posix()} does not exists.')

def getDescriptionText(config, relevant_configuration_keys):
    return '\n'.join([ f'{key} : {config[key]}' for key in relevant_configuration_keys ])


def concatenateVideos(
        paths,
        configs,
        relevant_configuration_keys,
        output_path,
        overwrite=False,
        font_color='black'):
    import ffmpeg

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

class Process(object):

    def __init__(self, output_directory, arguments):
        self.output_directory = output_directory
        self.arguments = arguments

    def getOutputDirectory(self, path):
        return self.output_directory / path

    def __call__(self, path_and_config_pairs):
        path, config_pairs = path_and_config_pairs

        experiment_path = self.getOutputDirectory(path)
        checkPathExists(experiment_path)

        config_path = experiment_path / 'config.txt'
        checkPathExists(config_path)
        config = SoapFilmSimulationConfigFile.fromConfigFile(config_path.as_posix())
        config.update(config_pairs)

        self.run(path, config)

        return (experiment_path, config)


class Render(Process):

    @staticmethod
    def getObjNumber(obj_path):
        return int(obj_path.name[4:10])

    def getFrameOutputPath(self, path):
        return pathlib.Path(self.arguments.frame_output_directory) / path / 'output'

    def getFrameFromObj(self, obj_path, path):
        return (
                self.getFrameOutputPath(path)
                / f'frame{Render.getObjNumber(obj_path)}.png'
            )

    @staticmethod
    def getObjs(path):
        output_path = path / 'output'
        for file in output_path.iterdir():
            if file.suffix == '.obj':
                yield file

    def renderObj(self, obj_path, path):
        output_frame = self.getFrameFromObj(obj_path, path)
        if self.arguments.skip_existing and output_frame.exists():
            return

        print(f'{obj_path} --> {output_frame}')
        if not self.arguments.test_run:
            subprocess.run(
                    [ 
                        self.arguments.blender_executable,
                        '--background',
                        '--python',
                        'render.py',
                        '--',
                        output_frame.as_posix(),
                        obj_path.as_posix(),
                        '--number-subdivisions', str(self.arguments.number_subdivisions_render),
                        '--scale', str(self.arguments.scale),
                        '--path-to-blender-toolbox', self.arguments.path_to_toolbox,
                        '--environment_texture', self.arguments.environment_texture
                    ],
                    capture_output=not self.arguments.show_blender_stdout
                )

    def run(self, path, config):
        experiment_path = self.getOutputDirectory(path)
        frame_output_path = self.getFrameOutputPath(path)
        if not frame_output_path.exists():
            frame_output_path.mkdir(parents=True)

        objs = list(Render.getObjs(experiment_path))
        objs.sort(key=Render.getObjNumber) 
        objs_and_path = [
                    (obj_path, path) for obj_path in objs[::self.arguments.render_every_n_obj]
                ]
        if self.arguments.number_threads > 1:
            list(pool.imap_unordered(
                self.renderObj,
                objs,
                max(len(objs) // (5 * self.arguments.number_threads), 1)))
            pool = multiprocessing.Pool(self.arguments.number_threads)
            list(pool.imap_unordered(
                self.renderObj,
                objs_and_path,
                max(len(objs) // (5 * self.arguments.number_threads), 1)))
            pool.close()
            pool.join()

        else:
            for obj, path in objs_and_path:
                self.renderObj(obj, path)
        


class Csv(Process):

    def run(self, path, config):
        experiment_path = self.getOutputDirectory(path)
        csv_path = experiment_path / 'data.csv'
        if csv_path.exists() and self.argument.skip_existing:
            return

        TOKEN = ['NumberVertices', 'FMMExecution', 'NaiveExecution']
        data = {}
        stdout = (experiment_path / 'stdout').read_text()
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
            csv_writer = csv.DictWriter(
                    csv_file,
                    fieldnames=data_keys,
                    delimiter=self.arguments.csv_delimiter)

            csv_writer.writeheader()
            csv_writer.writerows(data)

class Video(Process):

    def run(self, path, config):
        import ffmpeg
        experiment_path = self.getOutputDirectory(path)
        video_path = experiment_path / 'video.webm'
        if video_path.exists() and self.arguments.skip_existing:
            return

        
        temporary_frame_directory = experiment_path / 'Frames'
        divider = int(config.get('output-png-every-n-frames',
                    config.get('output-mesh-every-n-frames',
                        1)))

        if not self.arguments.skip_divide_frame_number:

            dividesFramesNumber(
                    'frame',
                    '.png',
                    temporary_frame_directory,
                    experiment_path / 'output',
                    divider)

        frame_rate = int(1 / (divider * config.get('time-step', 0.01)))

        stream = ffmpeg.input((temporary_frame_directory / 'frame%d.png').as_posix(), r=frame_rate)
        stream = ffmpeg.output(stream, video_path.as_posix())
        if self.arguments.overwrite:
            stream = ffmpeg.overwrite_output(stream)
        ffmpeg.run(stream)

class DeleteRender(Process):

    def getFrameOutputPath(self, path):
        return pathlib.Path(self.arguments.frame_output_directory) / 'output'

    def getFrameFromObj(self, obj_path, path):
        return (
                self.getFrameOutputPath(path)
                / f'frame{int(obj_path.name[4:10])}.png'
            )

    def run(self, path, config):
        for file in Render.getObjs(path):
            self.getFrameFromObj(file).unlink(missing_ok=True)


commands = ['video', 'render', 'delete-render', 'csv']

argument_parser = argparse.ArgumentParser()
argument_parser.add_argument('-o', '--sim-option', action='append', default=[])
argument_parser.add_argument('--directory-order', action='append', default=[])
argument_parser.add_argument('-j', '--number-threads', type=int, default=1)
argument_parser.add_argument('output_directory')
argument_parser.add_argument('command', choices=commands)
argument_parser.add_argument('command_arguments', nargs=argparse.REMAINDER)

commands_arguments_parser = {
        command : argparse.ArgumentParser(prog=sys.argv[0]+f' {command}')
        for command in commands 
        }

commands_arguments_parser['video'].add_argument(
        '--skip-existing', action='store_true')
commands_arguments_parser['video'].add_argument(
        '--overwrite', action='store_true')
commands_arguments_parser['video'].add_argument(
        '--font-color', default='black')
commands_arguments_parser['video'].add_argument(
        '--skip-divide-frame-number', action='store_true')

commands_arguments_parser['render'].add_argument(
        '--skip-existing', action='store_true')
commands_arguments_parser['render'].add_argument(
        '--scale', type=float, default=1)
commands_arguments_parser['render'].add_argument(
        '--path-to-toolbox', required=True)
commands_arguments_parser['render'].add_argument(
        '--number-subdivisions-render', type=int, default=0)
commands_arguments_parser['render'].add_argument(
        '--show-blender-stdout', action='store_true')
commands_arguments_parser['render'].add_argument(
        '--blender-executable', default='blender')
commands_arguments_parser['render'].add_argument(
        '-j', '--number_threads', type=int, default=1,
        help='This command number of threads should be kept to the default if the main program'
            'number of threads is different from the default')
commands_arguments_parser['render'].add_argument(
        '-t', '--test-run', action='store_true')
commands_arguments_parser['render'].add_argument(
        '--environment-texture')
commands_arguments_parser['render'].add_argument(
        '--frame-output-directory')
commands_arguments_parser['render'].add_argument(
        '-T', '--render-every-n-obj', type=int, default=1)

commands_arguments_parser['delete-render'].add_argument(
        '--frame-output-directory')

commands_arguments_parser['csv'].add_argument('--csv-delimiter', default=' ')
commands_arguments_parser['csv'].add_argument('--skip-existing', action='store_true')

args = argument_parser.parse_args()
command_args = commands_arguments_parser[args.command].parse_args(args.command_arguments)

output_directory = pathlib.Path(args.output_directory)
checkPathExists(output_directory)

simulation_parameter_product = SimulationParameterProduct(args.directory_order)
for sim_option in args.sim_option:
    key, values = sim_option.split('=')
    values = values.split(',')
    simulation_parameter_product.addOptions(key, values)

if args.command == 'video':
    process = Video(output_directory, command_args)
elif args.command == 'render':
    process = Render(output_directory, command_args)
elif args.command == 'delete-render':
    process = DeleteRender(output_directory, command_args)
elif args.command == 'csv':
    process = Csv(output_directory, command_args)

if args.number_threads > 1:
    with multiprocessing.Pool(args.number_threads) as pool:
        experiment_paths_and_config = list(pool.imap_unordered(
                process,
                [ (path, config_pairs) for path, config_pairs in simulation_parameter_product ]))
else:
    experiment_paths_and_config = []
    for path, config in simulation_parameter_product:
        experiment_paths_and_config.append(process((path, config)))

if args.command == 'video' and simulation_parameter_product.getNumberDifferentConfigurations() > 1:
    concatenateVideos(
            *zip(*experiment_paths_and_config),
            simulation_parameter_product.getRelevantConfigurationKeys(),
            output_directory / 'video.webm',
            font_color=command_args.font_color)

        
