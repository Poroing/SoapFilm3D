import argparse
import pathlib
import shutil

def dividesFramesNumber(
        prefix,
        suffix,
        output_dir,
        input_dir,
        divider,
        test_run=False,
        move=False,
        swap=False):

    input_directory_path = pathlib.Path(input_dir)
    if not input_directory_path.is_dir():
        raise ValueError(f'Input directory {input_dir.as_posix()} is not a directory')

    output_directory_path = pathlib.Path(output_dir)
    output_directory_path.mkdir(parents=True, exist_ok=True)

    for file_path in input_directory_path.iterdir():
        if not file_path.name.startswith(prefix) or not file_path.name.endswith(suffix):
            continue
        frame_number = int(file_path.name[len(prefix):-len(suffix)])
        if frame_number % divider != 0:
            continue
        new_frame_number = frame_number // divider
        
        if not swap:
            new_file_name = f'{prefix}{new_frame_number}{suffix}'
        else:
            new_file_name = f'{suffix}{new_frame_number}{prefix}'

        new_file_path = output_directory_path / new_file_name
        
        print(f'{file_path.as_posix()} --> {new_file_path.as_posix()}')
        if not test_run:
            if move:
                file_path.rename(new_file_path)
            else:
                shutil.copyfile(file_path.as_posix(), new_file_path.as_posix())
    


if __name__ == '__main__':
    argument_parser = argparse.ArgumentParser(
            description=
                'Divide the frame number of a set of numbered frame by a given number.\n'
                '\n'
                '   <SUFFIX>%d<PREFIX> -> <SUFFIX>(%d/n)<PREFIX>')

    argument_parser.add_argument('-p', '--prefix', default='frame')
    argument_parser.add_argument('-s', '--suffix', default='.png')
    argument_parser.add_argument('-o', '--output-dir', default='output')
    argument_parser.add_argument('-i', '--input-dir', default='.')
    argument_parser.add_argument('-t', '--test-run', action='store_true')
    argument_parser.add_argument('-m', '--move', action='store_true',
        help='Instead of making a copy of the file, moves the file')
    argument_parser.add_argument('--swap', action='store_true',
        help='This argument is here only to correct a mistake in the original version of the file.'
            'You don\'t need it.')
    argument_parser.add_argument('divider', type=int)
    args = argument_parser.parse_args()

    dividesFramesNumber(
            args.prefix,
            args.suffix,
            args.output_dir,
            args.input_dir,
            args.divider,
            test_run=args.test_run,
            move=args.move,
            swap=args.swap)
