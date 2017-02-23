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

    blast_hits = os.listdir(blast_dir)

    if len(blast_hits) < 1:
        print('Found no files in --blast_dir')
        exit(1)

    if not os.path.isfile(annot_db):
        print('--annot_db "{}" is not valid'.format(annot_db))
        exit(1)

    db = sqlite3.connect(annot_db)

    blast_fields = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch',
        'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
    sample_fields = ['sample_id', 'sample_acc', 'sample_name', 'project_name']
    sql = 'select ' + ', '.join(sample_fields) + ' from sample where sample_id=?'

    # print headers for output
    out_fh   = open(out_file, 'wt')
    out_fh.write('\t'.join(['qseqid', 'sample'] + sample_fields) + '\n')

    def err(msg):
        if args.verbose:
            sys.stderr.write(msg + '\n')

    for blast_out in blast_hits:
        match = re.match('.+-(\d+)\.tab$', blast_out)
        if not match:
            err('Failed to extract sample_id from blast-out "{}"'.format(blast_out))
            continue

        sample_id = match.group(1)
        for row in db.execute(sql, (sample_id,)):
            out_fh.write('\t'.join([blast_out] + list(map(str,row))) + '\n')

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
