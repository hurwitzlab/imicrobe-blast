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

That was too much. Instead build a sqlite3 database of
    FASTA file name
        sum of sequence lengths
        sum of seq_length * log(seq_length)
        sum of seq_length ** 2

This means parsing every FASTA file. Do that concurrently with 48 processes on Stampede2.

In addition this script writes a file of invalid FASTA files and a list of valid FASTA files.

test usage:
(imblst) jklynch@minty ~/host/project/imicrobe/apps/imicrobe-blast $ build_seq_db \
    -i "test/**/*.fa" \
    -d sqlite:///test_seq_db.sqlite \
    --invalid-files-fp test/bad_files.txt \
    --valid-files-fp test/good_files. \
    --work-dp test/work \
    --file-limit 1

"""
import argparse
import concurrent.futures
import glob
import itertools
import json
import os
import sys
import time

from Bio import SeqIO
from Bio.Alphabet import IUPAC

import numpy as np

from sqlalchemy import Column, Float, Integer, String, ForeignKeyConstraint, UniqueConstraint, create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship

from orminator import session_manager_from_db_uri


Base = declarative_base()


class FastaFile(Base):
    __tablename__ = 'fasta_file'

    id = Column(Integer, primary_key=True)

    file_path = Column(String, unique=True)

    seq_length_sum = Column(Float)
    seq_length_log_sum = Column(Float)
    seq_length_sq_sum = Column(Float)


class BadFastaFile(Base):
    __tablename__ = 'bad_fasta_file'

    id = Column(Integer, primary_key=True)

    file_path = Column(String, unique=True)
    exception_message = Column(String)


# class FastaSequence(Base):
#     __tablename__ = 'sequence'
#     __table_args__ = (ForeignKeyConstraint(('fasta_file_id',), ('fasta_file.id', )), )
#
#     id = Column(Integer, primary_key=True)
#
#     seq_id = Column(String)
#     seq_length = Column(Integer)
#     fasta_file_id = Column(Integer)
#
#     fasta_file = relationship('FastaFile')
#
#     UniqueConstraint('seq_id', 'fasta_file_id', name='unique_sequence_within_file')


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
    arg_parser.add_argument('--work-dp', default='work',
                            help='directory to write files')
    arg_parser.add_argument('--file-limit', type=int, required=False, default=None,
                            help='stop after <file-limit> files have been processed')

    args = arg_parser.parse_args(argv)
    print('command line arguments:\n\t{}'.format(args))

    return args


def main():
    build_seq_db(**vars(get_args(sys.argv[1:])))


def build_seq_db(fasta_globs, db_uri, invalid_files_fp, valid_files_fp, max_workers, work_dp, file_limit):
    """
    This method resulted in 'queue full' errors when run on all iMicrobe samples.

    Tried compressing the output of parse_fasta.

    :param fasta_globs:
    :param db_uri:
    :param invalid_files_fp:
    :param valid_files_fp:
    :param max_workers:
    :return:
    """
    fasta_list = []
    for fasta_glob in fasta_globs.split(','):
        glob_results = glob.glob(fasta_glob, recursive=True)
        print('glob "{}" matched {} files'.format(fasta_glob, len(glob_results)))
        fasta_list.extend((os.path.abspath(fp) for fp in glob_results))
    fasta_list.sort()

    print('SQLite db URL: {}'.format(db_uri))
    engine = create_engine(db_uri, echo=False)
    Base.metadata.create_all(engine, checkfirst=True)

    # get a list of all FASTA files already in the database
    # so they will not be loaded again
    with session_manager_from_db_uri(db_uri=db_uri) as db_session:
        fasta_fp_in_db = {fasta_file.file_path for fasta_file in db_session.query(FastaFile).all()}
        bad_fasta_fp_in_db = {bad_fasta_file.file_path for bad_fasta_file in db_session.query(BadFastaFile).all()}
        all_fasta_fp = set(fasta_list)
        all_fasta_fp_in_db = fasta_fp_in_db.union(bad_fasta_fp_in_db)
        fasta_fp_not_in_db = sorted(list(all_fasta_fp.difference(all_fasta_fp_in_db)))

    print('{} FASTA files found'.format(len(fasta_list)))
    print('{} FASTA file paths in database'.format(len(all_fasta_fp_in_db)))
    print('{} FASTA file paths not in database'.format(len(fasta_fp_not_in_db)))

    os.makedirs(work_dp, exist_ok=True)

    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        """Build a dict of future -> FASTA file path
        {
            future_1: '/work/05066/imicrobe/iplantc.org/data/ohana/HOT/HOT224_1_0025m/proteins.faa',
            future_2: '/work/05066/imicrobe/iplantc.org/data/ohana/HOT/HOT224_1_0050m/proteins.faa',
            ...
        }
        """
        future_to_fasta_fp = {
            executor.submit(parse_fasta, fasta_fp): fasta_fp
            for fasta_fp
            in itertools.islice(fasta_fp_not_in_db, file_limit)}

        for future in concurrent.futures.as_completed(future_to_fasta_fp):
            fasta_fp = future_to_fasta_fp[future]
            print('inserting sequence ids and lengths from {}'.format((fasta_fp)))
            try:
                seq_length_sum, seq_length_log_sum, seq_length_sq_sum, t = future.result()
                #with open(json_seq_id_to_seq_length_fp, 'rt') as f,\
                with session_manager_from_db_uri(db_uri=db_uri) as db_session:

                    #seq_id_to_seq_length = json.load(f)

                    t0 = time.time()
                    fasta_file = FastaFile(
                        file_path=fasta_fp,
                        seq_length_sum=seq_length_sum,
                        seq_length_log_sum=seq_length_log_sum,
                        seq_length_sq_sum=seq_length_sq_sum)
                    db_session.add(fasta_file)
                    #for n, (seq_id, seq_length) in enumerate(seq_id_to_seq_length.items()):
                    #    fasta_seq = FastaSequence(seq_id=seq_id, seq_length=seq_length)
                    #    fasta_seq.fasta_file = fasta_file
                    #    db_session.add(fasta_seq)

                    #    if (n + 1) % 1000 == 0:
                    #        db_session.flush()
                    #db_session.flush()
                    print('{:8.2f}s to insert sequence totals from "{}"'.format(
                        time.time()-t0,
                        fasta_fp))

            except FastaParseException as exc:
                with session_manager_from_db_uri(db_uri=db_uri) as db_session:
                    bad_fasta_file = BadFastaFile(file_path=fasta_fp, exception_message=str(exc))
                    db_session.add(bad_fasta_file)
                    #bad.append((fasta_fp, exc))
            #else:
            #    print('{:8.2f}s to parse "{}"'.format(t, fasta_fp))
            #    os.remove(json_seq_id_to_seq_length_fp)
            #    #good.append((fasta_fp, seq_id_to_seq_length, t))

    # log the numbers of valid and invalid files
    with session_manager_from_db_uri(db_uri=db_uri) as db_session:
        sorted_good = [
            fasta_file.file_path
            for fasta_file
            in db_session.query(FastaFile).order_by(FastaFile.file_path).all()]
        print('\n{} valid FASTA file(s)'.format(len(sorted_good)))
        with open(valid_files_fp, 'wt') as valid_file:
            valid_file.write('\n'.join(sorted_good))
            valid_file.write('\n')

        sorted_bad = db_session.query(BadFastaFile).order_by(BadFastaFile.file_path).all()
        print('{} invalid FASTA file(s)'.format(len(sorted_bad)))
        with open(invalid_files_fp, 'wt') as invalid_file:
            for bad_fasta_file in sorted_bad:
                invalid_file.write(bad_fasta_file.file_path)
                invalid_file.write('\n')
                sys.stdout.write(bad_fasta_file.file_path)
                sys.stdout.write('\n')
                sys.stdout.write(bad_fasta_file.exception_message)
                sys.stdout.write('\n')

    # what is in the db?
    with session_manager_from_db_uri(db_uri=db_uri) as db_session:
        file_count = db_session.query(FastaFile).count()
        print('inserted {} FASTA file(s)'.format(file_count))
        for f in db_session.query(FastaFile).all():
            print('FASTA file id: {}'.format(f.id))
            print('  file path: {}'.format(f.file_path))
            print('  sequence length sum         : {:5.2f}'.format(f.seq_length_sum))
            print('  sequence length log sum     : {:5.2f}'.format(f.seq_length_log_sum))
            print('  sequence length squared sum : {:5.2f}'.format(f.seq_length_sq_sum))

        #sequence_count = db_session.query(FastaSequence).count()
        #print('inserted {} sequence(s)'.format(sequence_count))
        #for s in db_session.query(FastaSequence).all():
        #    print('  sequence id: {}'.format(s.id))
        #    print('  sequence FASTA id: {}'.format(s.seq_id))
        #    print('  sequence length: {}'.format(s.seq_length))
        #    print('  sequence FASTA file id: {}'.format(s.fasta_file_id))


def parse_fasta(fasta_fp):
    t0 = time.time()

    seq_length_sum = 0.0
    seq_length_log_sum = 0.0
    seq_length_sq_sum = 0.0

    alphabet_dna = set(IUPAC.ambiguous_dna.letters)
    alphabet_protein = set(IUPAC.extended_protein.letters)
    # Ohana protein sequences often end in '*'
    alphabet_protein.add('*')
    seq_id_to_seq_length = {}
    for r, record in enumerate(SeqIO.parse(fasta_fp, format='fasta', alphabet=IUPAC.ambiguous_dna)):
        seq_letters = set(str(record.seq).upper())

        if len(record.seq) == 0:
            msg = '{}: Record {} has 0-length sequence\nid: {}'.format(fasta_fp, r+1, record.id)
            raise FastaParseException(msg)
        else:
            # all letters in record.seq must be in alphabet(s), but not all letters in alphabet(s) must in record.seq
            if len(seq_letters.difference(alphabet_dna)) == 0 or len(seq_letters.difference(alphabet_protein)) == 0:
                n = len(record.seq)

                seq_id_to_seq_length[record.id] = n
                seq_length_sum += n
                seq_length_log_sum += n * np.log(n)
                seq_length_sq_sum += n ** 2
            else:
                msg = '{}: Failed to parse sequence {}\nid: {}\nsequence: {}'.format(fasta_fp, r+1, record.id, record.seq[:1000])
                raise FastaParseException(msg)

    if len(seq_id_to_seq_length) == 0:
        msg = '{} is empty'.format(fasta_fp)
        raise FastaParseException(msg)

    t = time.time() - t0
    return seq_length_sum, seq_length_log_sum, seq_length_sq_sum, t


def write_parse_fasta(fasta_fp, work_dp):
    """
    Try to solve the MemoryError caused by too many worker results piling up in the queue.

    First compress the results. This helped a little.

    :param fasta_fp:
    :param worker_delay:
    :return:
    """

    json_fp = os.path.join(work_dp, os.path.basename(fasta_fp) + '.json')
    if os.path.exists(json_fp):
        try:
            t = 0
            with open(json_fp, 'rt') as fp:
                json.load(fp=fp)
            # fasta_fp has already been parsed
            print('file {} has already been parsed'.format(fasta_fp))
        except json.JSONDecodeError:
            # parsing did not finish
            print('rewriting file {}'.format(json_fp))
            seq_id_to_seq_length, t = parse_fasta(fasta_fp=fasta_fp)
            with open(json_fp, 'wt') as f:
                json.dump(seq_id_to_seq_length, f)
    else:
        print('parsing file {}'.format(fasta_fp))
        seq_id_to_seq_length, t = parse_fasta(fasta_fp=fasta_fp)
        with open(json_fp, 'wt') as f:
            json.dump(seq_id_to_seq_length, f)

    return json_fp, t


# def get_sequence_weights(db_uri):
#     """
#     Return a dictionary of file path to sequence read lengths:
#       {
#         '/path/to/file1.fasta': (100, 200, 50, ...),
#         '/path/to/file2.fasta': (150, 250, 100, ...),
#         ...
#       }
#     """
#
#     print('reading files and read lengths from {}'.format(db_uri))
#     t0 = time.time()
#     with session_manager_from_db_uri(db_uri=db_uri) as db_session:
#         file_path_to_read_lengths = {
#             f.file_path: query_sequence_weights(db_uri=db_uri, fasta_file_id=f.id)
#             for f
#             in db_session.query(FastaFile).all()}
#     print('done in {:5.2f}s'.format(time.time()-t0))
#
#     return file_path_to_read_lengths


# def query_sequence_weights(db_uri, fasta_file_id):
#     with session_manager_from_db_uri(db_uri=db_uri) as db_session:
#         return tuple([
#             s.seq_length
#             for s
#             in db_session.query(FastaSequence).filter(FastaSequence.fasta_file_id == fasta_file_id).all()])


# def get_sequence_weights_speedy(db_uri, max_workers):
#     """
#     :param db_uri:
#     :param max_workers:
#     :return:
#     """
#     with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
#         with session_manager_from_db_uri(db_uri=db_uri) as db_session:
#             """Build a dict of future -> FASTA file path
#             {
#                 future_1: '/work/05066/imicrobe/iplantc.org/data/ohana/HOT/HOT224_1_0025m/proteins.faa',
#                 future_2: '/work/05066/imicrobe/iplantc.org/data/ohana/HOT/HOT224_1_0050m/proteins.faa',
#                 ...
#             }
#             """
#             future_to_fasta_fp = {
#                     executor.submit(query_sequence_weights, db_uri, fasta_file.id): fasta_file.file_path
#                     for fasta_file
#                     in db_session.query(FastaFile).all()}
#
#             file_path_to_read_lengths = {
#                 future_to_fasta_fp[future]: future.result()
#                 for future
#                 in concurrent.futures.as_completed(future_to_fasta_fp)}
#
#     return file_path_to_read_lengths


if __name__ == '__main__':
    main()
