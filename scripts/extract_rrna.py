#!/usr/bin/env python3
"""Extract rRNA sequences from genomes based on RNAmmer result.

Usage:
    extract_rrna.py fnadir gffdir

Output:
    summary.tsv, 5s.fa, 16s.fa, 23s.fa
"""

__author__ = 'Qiyun Zhu'
__license__ = 'BSD-3-Clause'
__version__ = '0.0.1-dev'
__email__ = 'qiyunzhu@gmail.com'

from sys import argv, exit
from os import listdir

if len(argv) != 3:
    print("Usage: %s <sequence_dir> <RNAmmer_results_dir>")
    exit(1)

# genome sequence directory
fnadir = argv[1]

# RNAmmer result directory
gffdir = argv[2]

# rRNA types
genes = ('5s', '16s', '23s')

# reverse complement map
rc = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
      'Y': 'R', 'R': 'Y', 'S': 'S', 'W': 'W',
      'K': 'M', 'M': 'K', 'B': 'V', 'D': 'H',
      'H': 'D', 'V': 'B', 'N': 'N'}

# write header
with open('summary.tsv', 'w') as f:
    f.write('%s\n' % '\t'.join((
        'sample', 'operon', 'gene', 'score', 'start', 'end', 'strand')))

for file in listdir(gffdir):
    id = file.split('_', 1)[0]
    print(id)

    # read metadata
    hits = []
    with open('%s/%s' % (gffdir, file), 'r') as f:
        for line in f:
            line = line.rstrip('\r\n')
            if line.startswith('#'):
                continue
            # columns (GFF format): seqname, source, feature, start, end,
            # score, strand, frame, attribute
            x = line.split('\t')
            hits.append({'operon': x[0],
                         'start': int(x[3]),
                         'end': int(x[4]),
                         'score': float(x[5]),
                         'strand': x[6],
                         'gene': x[8].split('_')[0]})

    # read genome sequences
    nucls = set(x['operon'] for x in hits)
    nucls = {x: '' for x in nucls}
    cnucl = ''
    with open('%s/%s' % (fnadir, file.replace('.gff', '.fa')), 'r') as f:
        for line in f:
            line = line.rstrip('\r\n')
            if line.startswith('>'):
                x = line[1:].split()[0]
                cnucl = x if x in nucls else ''
            elif cnucl:
                nucls[cnucl] += line

    # output files
    fo = open('summary.tsv', 'a')
    fos = {x: open('%s.fa' % x, 'a') for x in genes}

    for hit in hits:

        # write summary
        fo.write('%s\t%s\t%s\t%.1f\t%d\t%d\t%s\n'
                 % (id, hit['operon'], hit['gene'], hit['score'],
                    hit['start'], hit['end'], hit['strand']))

        # write sequences
        seq = nucls[hit['operon']][hit['start'] - 1:hit['end']].upper()
        if hit['strand'] == '-':
            seq = ''.join([rc[nt] for nt in seq[::-1]])
        fos[hit['gene']].write('>%s_%s\n%s\n' % (id, hit['operon'], seq))

    fo.close()
    for gene in genes:
        fos[gene].close()
