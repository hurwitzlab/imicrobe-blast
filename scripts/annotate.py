#!/usr/bin/env python3
"""Annotate BLAST output"""

# Author: Ken Youens-Clark <kyclark@email.arizona.edu>

import argparse
import json
import os
import sqlite3

# --------------------------------------------------
def get_args():
    """get args"""
    parser = argparse.ArgumentParser(description='Annotate BLAST for iMicrobe')
    parser.add_argument('-b', '--blast_dir', help='BLAST out directory',
                        type=str, metavar='DIR', required=True)
    parser.add_argument('-a', '--annot_db', help='Annotation database',
                        type=str, metavar='DB',
                        default='/work/05066/imicrobe/iplantc.org/data/imicrobe-annotdb/annots.db')
    parser.add_argument('-o', '--out_file', help='Output file',
                        type=str, metavar='FILE',
                        default='blast-annotations.tab')
    return parser.parse_args()

# --------------------------------------------------
def main():
    """main"""
    args = get_args()
    out_file = args.out_file
    blast_dir = args.blast_dir
    annot_db = args.annot_db

    if not os.path.isdir(blast_dir):
        print('--blast_dir "{}" is not a directory'.format(blast_dir))
        exit(1)

    blast_hits = list(filter(lambda x: x.endswith('.tab'), os.listdir(blast_dir)))

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
    sql = 'select annots from annot where sample_id=?'
    annots = []
    fld_names = set()
    for blast_out in blast_hits:
        sample_id, _ = os.path.splitext(blast_out)
        for row in db.execute(sql, (sample_id,)):
            annot = json.loads(row[0])
            annots.append(annot)
            for key in annot.keys():
                fld_names.add(key)

    #
    # Print headers for output
    #
    fld_names.remove('sample_id')
    cols = ['sample_id'] + sorted(fld_names)
    out_fh = open(out_file, 'wt')
    out_fh.write('\t'.join(cols) + '\n')
    for annot in annots:
        vals = [annot.get(col) or 'NA' for col in cols]
        out_fh.write('\t'.join(vals) + '\n')

    out_fh.close()
    print('Done, see output file "{}"'.format(out_file))

# --------------------------------------------------
if __name__ == '__main__':
    main()
