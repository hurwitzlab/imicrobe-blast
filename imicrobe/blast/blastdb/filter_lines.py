"""
Read lines from standard input and pass them to standard output.
Filter out lines that match those in the file specified in the
first command line argument.
"""
import sys


def main():
    with open(sys.argv[1], 'rt') as file:
        filter_out_lines = {line for line in file}

    for line in sys.stdin:
        if line in filter_out_lines:
            filter_out_lines.remove(line)
        else:
            sys.stdout.write(line)


if __name__ == '__main__':
    main()