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

In addition this script writes a file of invalid FASTA files and a list of valid FASTA files.
"""
import argparse
import concurrent.futures
import glob
import os
import sys
import time

from Bio import SeqIO
from Bio.Alphabet import IUPAC

from sqlalchemy import Column, Integer, String, ForeignKeyConstraint, create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship

from orminator import session_manager_from_db_uri


Base = declarative_base()


class FastaFile(Base):
    __tablename__ = 'fasta_file'

    id = Column(Integer, primary_key=True)

    file_path = Column(String, unique=True)


class FastaSequence(Base):
    __tablename__ = 'sequence'
    __table_args__ = (ForeignKeyConstraint(('fasta_file_id',), ('fasta_file.id', )), )

    id = Column(Integer, primary_key=True)

    seq_id = Column(String, unique=True)
    seq_length = Column(Integer)
    fasta_file_id = Column(Integer)

    fasta_file = relationship('FastaFile')


class FastaParseException(BaseException):
    pass


def get_args(argv):
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument('-i', '--fasta-globs', required=True,
                            help='comma-separated list of globs for FASTA files to be processed')
    arg_parser.add_argument('-d', '--db-uri', required=True,
                            help='URI for sqlite3 database')
    arg_parser.add_argument('--invalid-files-fp', required=True,
                            help='path to output file of invalid FASTA files')
    arg_parser.add_argument('--valid-files-fp', required=True,
                            help='path to output file of valid FASTA files')
    arg_parser.add_argument('--max-workers', type=int, default=1,
                            help='number of processes')

    args = arg_parser.parse_args(argv)
    print('command line arguments:\n\t{}'.format(args))

    return args


def main():
    build_seq_db(**vars(get_args(sys.argv[1:])))


def build_seq_db(fasta_globs, db_uri, invalid_files_fp, valid_files_fp, max_workers):
    fasta_list = []
    for fasta_glob in fasta_globs.split(','):
        glob_results = glob.glob(fasta_glob, recursive=True)
        print('glob "{}" matched {} files'.format(fasta_glob, len(glob_results)))
        fasta_list.extend((os.path.abspath(fp) for fp in glob_results))

    print('SQLite db URL: {}'.format(db_uri))
    engine = create_engine(db_uri, echo=False)
    Base.metadata.create_all(engine)

    good = []
    bad = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        future_to_fasta_fp = {executor.submit(parse_fasta, fasta_fp): fasta_fp for fasta_fp in fasta_list}
        for future in concurrent.futures.as_completed(future_to_fasta_fp):
            fasta_fp = future_to_fasta_fp[future]
            try:
                seq_id_to_seq_length, t = future.result()
                with session_manager_from_db_uri(db_uri=db_uri) as db_session:
                    fasta_file = FastaFile(file_path=fasta_fp)
                    db_session.add(fasta_file)
                    for seq_id, seq_length in seq_id_to_seq_length.items():
                        fasta_seq = FastaSequence(seq_id=seq_id, seq_length=seq_length)
                        fasta_seq.fasta_file = fasta_file
                        db_session.add(fasta_seq)

            except FastaParseException as exc:
                bad.append((fasta_fp, exc))
            else:
                good.append((fasta_fp, seq_id_to_seq_length, t))

    sorted_good = sorted([fasta_fp for fasta_fp, *_ in good])
    print('\n{} valid FASTA file(s)\n'.format(len(sorted_good)))
    with open(valid_files_fp, 'wt') as valid_file:
        valid_file.write('\n'.join(sorted_good))
        valid_file.write('\n')

    sorted_bad = sorted([fasta_fp for (fasta_fp, _) in bad])
    print('{} invalid FASTA file(s):'.format(len(sorted_bad)))
    with open(invalid_files_fp, 'wt') as invalid_file:
        invalid_file.write('\n'.join(sorted_bad))
        invalid_file.write('\n')

    # what is in the db?
    with session_manager_from_db_uri(db_uri=db_uri) as db_session:
        for f in db_session.query(FastaFile).all():
            print('FASTA file id: {}'.format(f.id))
            print('FASTA file path: {}'.format(f.file_path))

        for s in db_session.query(FastaSequence).all():
            print('  sequence id: {}'.format(s.id))
            print('  sequence FASTA id: {}'.format(s.seq_id))
            print('  sequence length: {}'.format(s.seq_length))
            print('  sequence FASTA file id: {}'.format(s.fasta_file_id))


def parse_fasta(fasta_fp):
    t0 = time.time()
    alphabet = set(IUPAC.ambiguous_dna.letters)
    seq_id_to_seq_length = {}
    for r, record in enumerate(SeqIO.parse(fasta_fp, format='fasta', alphabet=IUPAC.ambiguous_dna)):
        seq_letters = set(str(record.seq).upper())

        if len(record.seq) == 0:
            msg = '{}: Record {} has 0-length sequence\nid: {}'.format(fasta_fp, r+1, record.id)
            raise FastaParseException(msg)
        else:
            # all letters in record.seq must be in alphabet, but not all letters in alphabet must in record.seq
            if len(seq_letters.difference(alphabet)) == 0:
                seq_id_to_seq_length[record.id] = len(record.seq)
            else:
                msg = '{}: Failed to parse sequence {}\nid: {}\nsequence: {}'.format(fasta_fp, r+1, record.id, record.seq[:1000])
                raise FastaParseException(msg)

    if len(seq_id_to_seq_length) == 0:
        msg = '{} is empty'.format(fasta_fp)
        raise FastaParseException(msg)

    t = time.time() - t0
    return seq_id_to_seq_length, t


if __name__ == '__main__':
    main()
