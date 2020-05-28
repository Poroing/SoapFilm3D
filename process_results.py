from foam_fmm_evaluation import *
import pathlib
import argparse
from divide_frame_number import dividesFramesNumber 
import csv
import sys
import multiprocessing

def hasTimedout(output_directory, path, set_time):
    return (output_directory / path / 'stdout').read_text().find(f'T = {float(set_time):g}') == -1

def hasError(output_directory, path):
    return (output_directory / path / 'error').exists()

def checkPathExists(path):
    if not path.exists():
        raise ValueError(f'{path} does not exists.')

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

def getStemInstance(config, parameters, stem):
    stem_instance = stem
    for key, value in config.items():
        if key in parameters:
            stem_instance = stem_instance + '-' + key + '=' + str(value)
    return stem_instance


def plot(simulation_parameter_product, args, output_directory):

    template_path = pathlib.Path(args.template)
    checkPathExists(template_path)
    template = template_path.read_text()

    for paths_and_configs in simulation_parameter_product.iterateOnSets(args.compare):

        output_filename = getStemInstance(
                paths_and_config[0][1],
                simulation_parameter_product.getRelevantConfigurationKeysNotInList(args.compare),
                args.output_file_prefix
            ) + '.tex'

        template_instance = template
        for path, config in paths_and_configs:

            csv_path = output_directory / path / (args.data_file_stem + '.csv')
            checkPathExists(csv_path)

            placeholder = '{' + getStemInstance(config, args.compare, args.placeholder) + '}'

            template_instance = template_instance.replace(placeholder, '{' + str(csv_path) + '}')

        output_file_path = output_directory / output_filename
        output_file_path.write_text(template_instance)

        subprocess.run([ 'latexmk', '-cd', '-pdf', output_file_path])

        if not args.no_convert_to_png:
            subprocess.run([
                    'convert',
                    '-density', '500',
                    output_file_path.with_suffix('.pdf'),
                    output_file_path.with_suffix('.png')
                ])

def compileExecutionTime(rows):
    execution_time = [ row.get('NaiveExecution', row.get('FMMExecution', None)) for row in rows ]
    number_vertices = [ row['NumberVertices'] for row in rows ]
    execution_time.sort()
    return {
            'MedianExecutionTime' : execution_time[len(execution_time) // 2],
            'MeanBiotSavartExecutionTime' : sum(map(float, execution_time)) / len(execution_time),
            'MeanNumberVertices' : sum(map(float, number_vertices)) / len(number_vertices)
        }

def compileMeanBubbleNumberVertices(rows):
    bubble_vertices = [ row['MeanBubbleNumberVertices'] for row in rows ]
    bubble_vertices.sort()
    return {
            "MedianBubbleVertices" : bubble_vertices[len(bubble_vertices) // 2],
            "MeanBubbleVertices" : sum(map(float, bubble_vertices)) / len(bubble_vertices)
        }

def compileStepExecutionTime(rows):
    step_execution_time = [ row['StepExecution'] for row in rows ]
    biot_savart_execution_time = [ row['BiotSavartExecution'] for row in rows ]
    number_vertices = [ row['NumberVertices'] for row in rows ]
    return {
            'MeanStepExecutionTime' : sum(map(float, step_execution_time)) / len(step_execution_time),
            'MeanStepBiotSavartExecutionTime' : sum(map(float, biot_savart_execution_time)) / len(biot_savart_execution_time),
        }

def compile(output_directory, path, stem, compile_function, args):
    data_path = output_directory / path / (stem + '.csv')
    with open(data_path, newline='') as data_file:
        data_object = csv.DictReader(data_file, delimiter=args.delimiter)
        compiled_data = compile_function(list(data_object))
    return compiled_data

def aggregate(
        simulation_parameter_product,
        args,
        output_directory,
        ignore_non_existing=False,
        ignore_timedout=False,
        ignore_error=False):

    
    for paths_and_configs in simulation_parameter_product.iterateOnSets([ args.aggregate_on ]):
        data = []
        for path, config in paths_and_configs:
            if not (output_directory / path).exists() and ignore_non_existing:
                continue
            if (
                    hasTimedout(output_directory, path, config.get('simulation-time', 4.0))
                    and ignore_timedout
                ):
                continue
            if hasError(output_directory, path) and ignore_error:
                continue

            compiled_data = { args.aggregate_on : config[args.aggregate_on] }
            if args.execution_time_file_stem is not None:
                compiled_data.update(compile(
                        output_directory,
                        path,
                        args.execution_time_file_stem,
                        compileExecutionTime,
                        args
                    ))
            if args.mean_number_vertices_file_stem is not None:
                compiled_data.update(compile(
                        output_directory,
                        path,
                        args.mean_number_vertices_file_stem,
                        compileMeanBubbleNumberVertices,
                        args
                    ))
            if args.step_execution_time_file_stem is not None:
                compiled_data.update(compile(
                        output_directory,
                        path,
                        args.step_execution_time_file_stem,
                        compileStepExecutionTime,
                        args
                    ))
            data.append(compiled_data)

        output_filename = getStemInstance(
                paths_and_configs[0][1],
                simulation_parameter_product.getRelevantConfigurationKeysNotInList(
                        [ args.aggregate_on ]
                    ),
                args.output_file_stem
            ) + '.csv'
        with open(output_directory / output_filename, newline='', mode='w') as output_file:
            output_object = csv.DictWriter(output_file, data[0].keys(), delimiter=args.delimiter)
            output_object.writeheader()
            output_object.writerows(data)


class Process(object):

    def __init__(self, output_directory, arguments):
        self.output_directory = output_directory
        self.arguments = arguments
        self.ignore_non_existing = False

    def getOutputDirectory(self, path):
        return self.output_directory / path

    def __call__(self, path_and_config_pairs):
        path, config_pairs = path_and_config_pairs

        experiment_path = self.getOutputDirectory(path)
        if self.ignore_non_existing and not experiment_path.exists():
            return (None, None)

        checkPathExists(experiment_path)

        config_path = experiment_path / 'config.txt'
        if config_path.exists():
            config = SoapFilmSimulationConfigFile.fromConfigFile(config_path.as_posix())
        else:
            config = SoapFilmSimulationConfigFile()
        config.update(config_pairs)

        self.run(path, config)

        return (experiment_path, config)


class Render(Process):

    @staticmethod
    def getObjNumber(obj_path):
        return int(obj_path.name[4:10])

    def getFrameOutputPath(self, path):
        if self.arguments.frame_output_directory is None:
            return self.getOutputDirectory(path) / 'output'
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

    def renderObj(self, obj_path_and_path):
        obj_path, path = obj_path_and_path
        output_frame = self.getFrameFromObj(obj_path, path)
        if self.arguments.skip_existing and output_frame.exists():
            print(f'{output_frame} exists, skipping.')
            return

        print(f'{obj_path} --> {output_frame}')
        if not self.arguments.test_run:
            command = [ 
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
                        ]
            
            if self.arguments.environment_texture is not None:
                command.extend([
                        '--environment_texture', self.arguments.environment_texture
                        ])

            subprocess.run(
                    command,
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
            pool = multiprocessing.Pool(self.arguments.number_threads)
            list(pool.imap_unordered(
                self.renderObj,
                objs_and_path,
                max(len(objs) // (5 * self.arguments.number_threads), 1)))
            pool.close()
            pool.join()

        else:
            for obj_and_path in objs_and_path:
                self.renderObj(obj_and_path)
        
class Csv(Process):

    def run(self, path, config):
        experiment_path = self.getOutputDirectory(path)
        csv_path = experiment_path / (self.arguments.filename_stem + '.csv')
        if csv_path.exists() and self.arguments.skip_existing:
            return

        data = {}
        stdout = (experiment_path / 'stdout').read_text()
        for line in stdout.split('\n'):
            line_tokens = line.split(' ')
            line_data = self.processLine(line_tokens)
            if line_data is not None:
                for key, value in line_data.items():
                    data.setdefault(key, []).append(value)

        data_length = None
        for value in data.values():
            if data_length is None:
                data_length = len(value)

            if data_length != len(value) and not self.arguments.crop:
                raise ValueError('The standard output does not contain as much data of each type')

            data_length = min(len(value), data_length)

        data_keys = data.keys()
        data = [ { k : v[i] for k, v in data.items() } for i in range(data_length) ]
        with open(csv_path.as_posix(), mode='w', newline='') as csv_file:
            csv_writer = csv.DictWriter(
                    csv_file,
                    fieldnames=data_keys,
                    delimiter=self.arguments.csv_delimiter)

            csv_writer.writeheader()
            csv_writer.writerows(data)

class CsvExecutionTime(Csv):
    TOKEN = ['NumberVertices', 'FMMExecution', 'NaiveExecution']

    def processLine(self, line_tokens):
        if line_tokens[0] not in self.TOKEN:
            return None

        return { line_tokens[0] : float(line_tokens[1]) }

class CsvMeanNumberVertices(Csv):

    def processLine(self, line_tokens):
        if line_tokens[0] != 'NumberVerticesIncidentToRegions':
            return None

        bubble_number_vertices = list(map(int, line_tokens[2:]))
        mean_bubble_number_vertices = sum(bubble_number_vertices) / len(bubble_number_vertices)
        return { 'MeanBubbleNumberVertices' : mean_bubble_number_vertices }

class CsvTimeStepExecutionTime(CsvExecutionTime):

    def run(self, path, config):
        self.total_biot_savart_execution_time = 0
        self.number_vertices = 0
        super().run(path, config)


    def processLine(self, line_tokens):
        data = super().processLine(line_tokens)
        if data is not None:
            if 'NumberVertices' in data:
                self.number_vertices = data['NumberVertices']
            elif 'NaiveExecution' in data:
                self.total_biot_savart_execution_time += data['NaiveExecution']
            elif 'FMMExecution' in data:
                self.total_biot_savart_execution_time += data['FMMExecution']
            return None
        elif line_tokens[0] == 'StepExecution':
            data = {
                'NumberVertices' : self.number_vertices,
                'BiotSavartExecution' : self.total_biot_savart_execution_time,
                'StepExecution' : float(line_tokens[1])
            }
            self.total_biot_savart_execution_time = 0
            return data
        return None






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

commands = ['video', 'render', 'delete-render', 'csv', 'plot', 'aggregate']

argument_parser = argparse.ArgumentParser()
argument_parser.add_argument('-o', '--sim-option', action='append', default=[])
argument_parser.add_argument('--directory-order', action='append', default=[])
argument_parser.add_argument('-j', '--number-threads', type=int, default=1)
argument_parser.add_argument('--ignore-non-existing', action='store_true')
argument_parser.add_argument('--ignore-timeout', action='store_true')
argument_parser.add_argument('--ignore-error', action='store_true')
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
commands_arguments_parser['video'].add_argument(
        '--no-concatenation', action='store_true')

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

commands_arguments_parser['csv'].add_argument(
        '--csv-delimiter', default=' ')
commands_arguments_parser['csv'].add_argument(
        '--skip-existing', action='store_true')
commands_arguments_parser['csv'].add_argument(
        '--mode',
        choices=['mean_number_vertices', 'execution_time', 'step_execution_time'],
        required=True)
commands_arguments_parser['csv'].add_argument(
        '--filename-stem', default='data')
commands_arguments_parser['csv'].add_argument(
        '--crop', action='store_true')

commands_arguments_parser['plot'].add_argument(
        '--template', required=True)
commands_arguments_parser['plot'].add_argument(
        '--placeholder', default='placeholder')
commands_arguments_parser['plot'].add_argument(
        '--no-convert-to-png', action='store_true')
commands_arguments_parser['plot'].add_argument(
        '--compare', action='append', default=[])
commands_arguments_parser['plot'].add_argument(
        '--output-file-prefix', default='plot')
commands_arguments_parser['plot'].add_argument(
        '--data-file-stem', default='data')

commands_arguments_parser['aggregate'].add_argument(
        '--aggregate-on', required=True)
commands_arguments_parser['aggregate'].add_argument(
        '--execution-time-file-stem')
commands_arguments_parser['aggregate'].add_argument(
        '--mean-number-vertices-file-stem')
commands_arguments_parser['aggregate'].add_argument(
        '--delimiter', default=' ')
commands_arguments_parser['aggregate'].add_argument(
        '--output-file-stem', default='aggregate')
commands_arguments_parser['aggregate'].add_argument(
        '--step-execution-time-file-stem')

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
elif args.command == 'csv' and command_args.mode == 'mean_number_vertices':
    process = CsvMeanNumberVertices(output_directory, command_args)
elif args.command == 'csv' and command_args.mode == 'execution_time':
    process = CsvExecutionTime(output_directory, command_args)
elif args.command == 'csv' and command_args.mode == 'step_execution_time':
    process = CsvTimeStepExecutionTime(output_directory, command_args)
else:
    process = None

if process is not None:
    to_process = []
    for path, config in simulation_parameter_product:
        if not (output_directory / path).exists() and args.ignore_non_existing:
            continue
        if (
                hasTimedout(output_directory, path, config.get('simulation-time', 4.0))
                and args.ignore_timeout
            ):
            continue
        if hasError(output_directory, path) and args.ignore_error:
            continue

        to_process.append((path, config))

    if args.number_threads > 1:
        with multiprocessing.Pool(args.number_threads) as pool:
            experiment_paths_and_config = list(pool.imap_unordered(
                    process,
                    to_process))

    else:
        experiment_paths_and_config = []
        for path, config in to_process:
            experiment_paths_and_config.append(process((path, config)))

if (
        args.command == 'video'
        and simulation_parameter_product.getNumberDifferentConfigurations() > 1
        and not command_args.no_concatenation
    ):
    concatenateVideos(
            *zip(*experiment_paths_and_config),
            simulation_parameter_product.getRelevantConfigurationKeys(),
            output_directory / 'video.webm',
            font_color=command_args.font_color)

if args.command == 'plot':
    plot(simulation_parameter_product, command_args, output_directory)

if args.command == 'aggregate':
    aggregate(simulation_parameter_product, command_args, output_directory,
            ignore_error=args.ignore_error, ignore_non_existing=args.ignore_non_existing,
            ignore_timedout=args.ignore_timeout)

        
