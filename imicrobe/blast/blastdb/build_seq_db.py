"""
Build a sqlite3 database of
    FASTA file name
        sequence id
        sequence length
        sequence id
        sequence length
        ...
        sequence id
        sequence length

This means parsing every FASTA file. Do that concurrently with 48 processes on Stampede2.
"""
import argparse
import concurrent.futures
import glob
import sys
import time

from Bio import SeqIO
from Bio import Alphabet
from Bio.Alphabet import IUPAC

from sqlalchemy import Column, String
from sqlalchemy import create_engine

from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker

from ormniator import session_from_uri


Base = declarative_base()


class FastaFile(Base):
    __tablename__ = 'fasta_file'

    id = Column(Integer, primary_key=True)

    file_path = Column(String)


class FastaSequence(Base):
    __tablename__ = 'sequence'

    id = Column(String, primary_key=True)

    seq_length = Column(Integer)

    fasta_file_id = Column()


def get_args(argv):
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument('-i', '--fasta-glob', required=True, help='glob for FASTA files to be validated')
    arg_parser.add_argument('-d', '--db', required=true, help='file path for sqlite3 database')
    arg_parser.add_argument('--max-workers', type=int, default=1, help='number of processes')

    args = arg_parser.parse_args(argv)
    print('command line arguments:\n\t{}'.format(args))

    return args


def main():
    build_seq_db(**vars(get_args(sys.argv[1:])))


def build_seq_db(fasta_glob, db, max_workers):
    fasta_list = glob.glob(fasta_glob, recursive=True)
    print('glob "{}" matched {} files'.format(fasta_glob, len(fasta_list)))

    good = []
    bad = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        future_to_fasta_fp = {executor.submit(parse_fasta, fasta_fp): fasta_fp for fasta_fp in fasta_list}
        for future in concurrent.futures.as_completed(future_to_fasta_fp):
            fasta_fp = future_to_fasta_fp[future]
            try:
                seq_id_to_seq_length, t = future.result()
                with
            except Exception as exc:
                bad.append((fasta_fp, exc))
            else:
                good.append((fasta_fp, read_count, t))

        print('\n{} valid FASTA file(s)\n'.format(len(good)))

        sorted_bad = sorted(bad)
        print('{} problematic file(s):'.format(len(bad)))
        print('\n'.join([b[0] for b in sorted_bad]))

        print('\nfailures:')
        print('\n\n'.join([str(b[1]) for b in sorted_bad]))


def parse_fasta(fasta_fp):
    t0 = time.time()
    alphabet = set(IUPAC.ambiguous_dna.letters)
    seq_id_to_seq_length = {}
    for record in SeqIO.parse(fasta_fp, format='fasta', alphabet=IUPAC.ambiguous_dna):
        seq_letters = set(str(record.seq).upper())

        if len(record.seq) == 0:
            msg = '{}: Record {} has 0-length sequence\nid: {}'.format(fasta_fp, read_count, record.id)
            raise Exception(msg)
        else:
            # all letters in record.seq must be in alphabet, but not all letters in alphabet must in record.seq
            if len(seq_letters.difference(alphabet)) == 0:
                seq_id_to_seq_length[record.id] = len(record.seq)
            else:
                msg ='{}: Failed to parse sequence {}\nid: {}\nsequence: {}'.format(fasta_fp, read_count, record.id, record.seq[:1000])
                raise Exception(msg)

    if len(seq_id_to_seq_length) == 0:
        msg = '{} is empty'.format(fasta_fp)
        raise Exception(msg)

    t = time.time() - t0
    return seq_id_to_seq_length, t


if __name__ == '__main__':
    main()
