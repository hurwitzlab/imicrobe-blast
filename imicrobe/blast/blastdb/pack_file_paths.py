"""
"""
import argparse
import itertools
import operator
import os.path
import pprint
import sys

import numpy as np


def get_args():
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument('-n', '--split-count', type=int, required=True,
                            help='Number of files into which stdin will be split')
    arg_parser.add_argument('--prefix', required=True,
                            help='File path prefix for split files')

    return arg_parser.parse_args()


def pack_file_lists(file_list, bin_count):
    """
    Given 'file_list' of (file size, file path) tuples return a list of 'bin_count' lists
    of file paths such that the files in each list have approximately the same total number
    of bytes.

    This function uses the "first-fit decreasing algorithm" described on Wikipedia's
    "bin packing problem" page: https://en.wikipedia.org/wiki/Bin_packing_problem.

    Arguments:
    file_list: sequence of tuples (file size, file path)
    bin_count: number of bins into which files will be packed
    """
    target_bin_weight = np.round(np.sum([weight for weight, name in file_list]) / bin_count)
    print('the target bin weight is {:8.1f}'.format(target_bin_weight))

    bin_list = [dict((('bin_weight', 0.0), ('bin_contents', []))) for _ in range(bin_count)]

    for (weight, name) in sorted(file_list, reverse=True):
        bin_list[0]['bin_contents'].append((weight, name))
        bin_list[0]['bin_weight'] += weight
        bin_list.sort(key=lambda bin_: bin_['bin_weight'])

    print('minimum bin weight: {:8.1f}\nmaximum bin weight: {:8.1f}'.format(
        bin_list[0]['bin_weight'], bin_list[-1]['bin_weight']))

    return [bin_['bin_contents'] for bin_ in bin_list]


def make_packed_file_lists(file_paths, file_list_count):
    """
    """
    print('{} file paths will be packed into {} lists'.format(len(file_paths), file_list_count))

    file_size_path_list = [(os.path.getsize(fp), fp) for fp in file_paths]
    print('first 5 files:\n{}'.format(pprint.pformat(file_size_path_list[:5])))

    packed_file_lists = pack_file_lists(file_size_path_list, bin_count=file_list_count)

    # remove the file weights and sort by file path
    sorted_packed_file_lists = [
        sorted([fp for weight, fp in packed_file_list])
        for packed_file_list
        in packed_file_lists
    ]

    return sorted_packed_file_lists


def main():
    args = get_args()

    if args.split_count < 2:
        print('--split-count must be greater than 1')
        quit()

    file_paths = sys.stdin.readlines()
    packed_file_lists = make_packed_file_lists(file_paths=file_paths, file_list_count=args.split_count)

    # this iterator yields 'aa', 'ab', 'ac', ..., 'zz'
    group_id_iter = itertools.starmap(
        operator.add,
        # this iterator yields ['a', 'a'], ['a', 'b'], ['a', 'c'], ..., ['z', 'z']
        itertools.product('abcdefghijklmnopqrstuvwxyz', repeat=2))

    for group_id, packed_file_list in zip(group_id_iter, packed_file_lists):
        file_name = args.prefix + group_id
        with open(file_name, 'w') as split_file:
            split_file.write(''.join(packed_file_list))


if __name__ == '__main__':
    main()
