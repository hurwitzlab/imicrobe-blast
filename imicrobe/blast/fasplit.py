# Author: Ken Youens-Clark <kyclark@email.arizona.edu>
# Author: Joshua Lynch <jklynch@email.arizona.edu>

import argparse
import itertools
import os

from Bio import SeqIO


def main():
    args    = get_args()
    fasta   = args.fasta
    out_dir = args.out_dir
    nfile   = args.splits

    if not os.path.isfile(fasta):
        print('--fasta "{}" is not valid'.format(fasta))
        exit(1)

    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    if nfile < 1:
        print("--splits cannot be less than one")
        exit(1)

    nseq = 0
    basename, ext = os.path.splitext(os.path.basename(fasta))

    try:
        # open all output files
        files = [open(os.path.join(out_dir, '{}_{}{}'.format(basename, file_number, ext), 'wt')) for file_number in range(nfile)]
        files_cycle = itertools.cycle(files)

        for record in SeqIO.parse(fasta, "fasta"):
            SeqIO.write(record, files_cycle.next(), "fasta")
            nseq += 1
    finally:
        for f in files:
            f.close()

    print('Done, wrote {} sequence{} to {} file{}'.format(
        nseq, '' if nseq == 1 else 's',
        nfile, '' if nfile == 1 else 's'))


def get_args():
    parser = argparse.ArgumentParser(description='Split FASTA files')
    parser.add_argument('-f', '--fasta', help='FASTA input file',
        type=str, metavar='FILE', required=True)
    parser.add_argument('-n', '--splits', help="Split the input into this many files",
        type=int, required=True)
    parser.add_argument('-o', '--out-dir', help='Output directory',
        type=str, metavar='DIR', default='fasplit')
    return parser.parse_args()


if __name__ == '__main__':
    main()
