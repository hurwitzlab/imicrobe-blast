"""
Testing:

(imblst) jklynch@minty ~/host/project/imicrobe/apps/imicrobe-blast $ time pack_file_paths -d sqlite:///ohana_little_seq_db.sqlite -n 2 --prefix bbb --max-workers 2

"""
import argparse
import itertools
import operator
import pprint
import time

import numpy as np

from imicrobe.blast.blastdb.build_seq_db import get_sequence_weights_speedy


def get_args():
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument('-n', '--split-count', type=int, required=True,
                            help='Number of files into which stdin will be split')
    arg_parser.add_argument('-d', '--db-uri', required=True,
                            help='URI for sqlite3 database')
    arg_parser.add_argument('--prefix', required=True,
                            help='File path prefix for split files')
    arg_parser.add_argument('--max-workers', type=int, default=1,
                            help='number of processes')

    return arg_parser.parse_args()


def pack_file_lists(file_list, bin_count):
    """
    Given 'file_list' of (file weight, file path) tuples return a list of 'bin_count' lists
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

    bin_list = [dict((('weight', 0.0), ('contents', []))) for _ in range(bin_count)]

    for (weight, name) in sorted(file_list, reverse=True):
        bin_list[0]['contents'].append((weight, name))
        bin_list[0]['weight'] += weight
        bin_list.sort(key=lambda bin_: bin_['weight'])

    print('minimum bin weight: {:8.1f}\nmaximum bin weight: {:8.1f}'.format(
        bin_list[0]['weight'], bin_list[-1]['weight']))

    return bin_list


def make_packed_file_lists(file_size_path_list, file_list_count):
    """Return a list of dictionaries:
        [
            {
                'contents': [('/file/path/1', 123456), ('/file/path/2', 234567), ...],
                'weight': 2000000
            },
            {
                'contents': [('/file/path/3', 345678), ('/file/path/4', 45678), ...],
                'weight': 3000000
            },
            ...
        ]
    """
    print('{} file paths will be packed into {} lists'.format(
        len(file_size_path_list), file_list_count))

    print('first 5 files:\n{}'.format(pprint.pformat(file_size_path_list[:5])))
    print('file weights:\n{}'.format(pprint.pformat(file_size_path_list)))

    packed_file_lists = pack_file_lists(file_size_path_list, bin_count=file_list_count)
    for packed_files in packed_file_lists:
        print(packed_files)
        print('weight: {:8.1f}\n\t{}'.format(
            packed_files['weight'],
            '\n\t'.join(['{:8.1f} {}'.format(weight, name) for (weight, name) in packed_files['contents']])))

    return packed_file_lists


def main():
    args = get_args()

    if args.split_count < 2:
        print('--split-count must be greater than 1')
        quit()

    t0 = time.time()
    #file_paths_to_sequence_lengths = get_sequence_weights(args.db_uri)
    #print('{:5.2f}s for get_sequence_weights'.format(time.time()-t0))
    file_paths_to_sequence_lengths = get_sequence_weights_speedy(args.db_uri, max_workers=args.max_workers)
    print('{:5.2f}s for get_sequence_weights_speedy'.format(time.time()-t0))

    t0 = time.time()
    packed_file_lists = make_packed_file_lists(
        file_size_path_list=tuple([
            (np.sum([(n**2)/1000 for n in read_length_list]), file_path)
            for file_path, read_length_list
            in file_paths_to_sequence_lengths.items()
        ]),
        file_list_count=args.split_count)
    print('{:5.2f}s for make_packed_file_lists done'.format(time.time()-t0))
    print('packed lists:')

    # this iterator yields 'aa', 'ab', 'ac', ..., 'zz'
    group_id_iter = itertools.starmap(
        operator.add,
        # this iterator yields ['a', 'a'], ['a', 'b'], ['a', 'c'], ..., ['z', 'z']
        itertools.product('abcdefghijklmnopqrstuvwxyz', repeat=2))

    for group_id, packed_file_list in zip(group_id_iter, packed_file_lists):
        file_name = args.prefix + group_id
        with open(file_name, 'w') as split_file:
            split_file.write('\n'.join([file_path for (file_path, weight) in packed_file_list['contents']]))


if __name__ == '__main__':
    main()
