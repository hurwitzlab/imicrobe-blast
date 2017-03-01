#!/usr/bin/env python

# Author: Ken Youens-Clark <kyclark@email.arizona.edu>

import argparse
import os
import re
import sqlite3
import sys

def main():
    args      = get_args()
    out_file  = args.out_file
    blast_dir = args.blast_dir
    annot_db  = args.annot_db
    verbose   = args.verbose

    if not os.path.isdir(blast_dir):
        print('--blast_dir "{}" is not a directory'.format(blast_dir))
        exit(1)

    blast_hits = filter(lambda x: x.endswith('.tab'), os.listdir(blast_dir))

    if len(blast_hits) < 1:
        print('Found no files in --blast_dir')
        exit(1)

    if not os.path.isfile(annot_db):
        print('--annot_db "{}" is not valid'.format(annot_db))
        exit(1)

    out_dir = os.path.dirname(out_file)
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    db = sqlite3.connect(annot_db)

    #
    # Find the field names (2nd in tuple returned by pragma)
    #
    cursor = db.execute('pragma table_info(sample)');
    sample_fields = map(lambda x: x[1], cursor.fetchall()) 
    sql = 'select * from sample where sample_id=?'

    #
    # Print headers for output
    #
    out_fh = open(out_file, 'wt')
    out_fh.write('\t'.join(sample_fields) + '\n')

    for blast_out in blast_hits:
        sample_id, ext = os.path.splitext(blast_out)
        for row in db.execute(sql, (sample_id,)):
            out_fh.write('\t'.join(list(map(str,row))) + '\n')

    out_fh.close()
    print('Done, see output file "{}"'.format(out_file))

def get_args():
    parser = argparse.ArgumentParser(description='Annotate BLAST for iMicrobe')
    parser.add_argument('-b', '--blast_dir', help='BLAST out directory',
        type=str, metavar='DIR', required=True)
    parser.add_argument('-a', '--annot_db', help='Annotation database',
        type=str, metavar='FILE',
        default='/work/03137/kyclark/imicrobe/annotations/')
    parser.add_argument('-o', '--out_file', help='Output file',
        type=str, metavar='DIR', default='blast-annotated')
    parser.add_argument('-v', '--verbose', help='Say more stuff',
        action='store_true')
    return parser.parse_args()

if __name__ == '__main__':
    main()
