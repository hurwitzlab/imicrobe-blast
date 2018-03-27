"""
Read all sequences in all FASTA files specified on standard input.
For each file write a catalog of the sequences that looks like this:

file_name.fa.catalog:
    id    length
    abc   53
    def   149
    ...   ...

Then calculate the weight for each file and pass weight and file path on
to standard output.
"""
import glob
import sys
import time

from Bio import SeqIO


with open('cache.txt', 'wt') as cache_file:

    for fasta_fp in (line.strip() for line in sys.stdin.readlines()):
        cache_file.write(fasta_fp)
        cache_file.write('\n')

        #_, fasta_filename = os.path.split(fasta_fp)

        weight = 0.0
        t0 = time.time()
        for record in SeqIO.parse(fasta_fp, format='fasta'):
            n = len(record.seq)
            weight += n ** 2
            cache_file.write(str(n))
            cache_file.write('\t')
            cache_file.write(record.id)
            cache_file.write('\n')
        cache_file.write('\n')

        print('parsed file {} in {:5.2f}s'.format(fasta_fp, time.time()-t0))
        print('  weight is {}'.format(round(weight)))

        sys.stdout.write('{},{}\n'.format(round(weight), fp))
