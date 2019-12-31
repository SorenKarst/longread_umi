#!/usr/bin/env python3
"""Assign taxonomy based on BLAST result.
"""

__author__ = 'Qiyun Zhu'
__license__ = 'BSD-3-Clause'
__version__ = '0.0.1-dev'
__email__ = 'qiyunzhu@gmail.com'

import sys
import argparse

usage = """%(prog)s -i INPUT_TABLE -t TAXONOMY -o OUTPUT_MAP [options]"""

description = """example:
  %(prog)s -i blast.out -t reftax.txt -k 0.1 -o taxonomy.txt
"""

epilog = """This script assigns taxonomy to each query sequence based on the
consensus taxonomy of top hit(s). It assumes that queries appear sequentially
and hits per query are ordered by bit score from high to low.
"""


def parse_args():
    parser = argparse.ArgumentParser(
        usage=usage, description=description, epilog=epilog,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    arg = parser.add_argument
    arg('-i', '--input', type=argparse.FileType('r'), default=sys.stdin,
        help='input BLAST hit table in tabular format (-outfmt 6)')
    arg('-o', '--output', type=argparse.FileType('w'), default=sys.stdout,
        help='output taxonomy map, format: query <tab> taxonomic unit')
    arg('-t', '--taxonomy', type=argparse.FileType('r'), required=True,
        help=('reference taxonomy, format: subject <tab> taxonomic;units;'
              'from;high;to;low"'))
    arg('-k', '--top', type=float, default=1,
        help=('top hits to consider; if k > 1, take top k hits; if 0 < k < 1,'
              ' select hits with bit score within k fraction below top hit.'))
    for arg in parser._actions:
        arg.metavar = '\b'
    return parser.parse_args()


def main():
    # parser arguments
    args = parse_args()

    # validate parameter
    if args.top > 1:
        args.top = int(args.top)
    elif args.top <= 0:
        raise ValueError('Invalid top hits cutoff: %s.' % args.top)

    # read reference taxonomy
    reftax = dict(x.split('\t') for x in args.taxonomy.read().splitlines())

    cqry = None  # current query
    csubs = []   # current subjects
    cbit = None  # bit score of top hit of current query
    skip = False  # skip current hit

    # assign taxonomy based on top hits
    def _assign_tax():
        nonlocal cqry, csubs, cbit, skip
        if cqry is None:
            return

        # simply adopt taxonomy if only one hit
        elif len(csubs) == 1:
            args.output.write('%s\t%s\n' % (cqry, reftax[csubs[0]]))

        # find consensus taxonomy of multiple hits
        else:
            clin = None  # current lineage
            for taxon in reftax[csubs[0]].split('; '):
                # take top hit as reference
                lin_ = taxon if not clin else '%s; %s' % (clin, taxon)

                # check all remaining hits for consistency
                npass = 1
                for sub in csubs[1:]:
                    tax = reftax[sub]
                    if tax == lin_ or tax.startswith('%s; ' % lin_):
                        npass += 1

                # halt if inconsistency is observed
                if npass < len(csubs):
                    break

                # otherwise extend lineage
                clin = lin_

            # output assignment
            args.output.write('%s\t%s\n' % (cqry, clin or 'Unassigned'))

        # clear intermediates
        csubs, skip = [], True

    # process BLAST hit table
    for line in args.input:
        x = line.rstrip('\r\n').split('\t')
        if x[0] == cqry:
            if skip:
                continue

            # stop when k top hits have been selected
            if args.top >= 1 and len(csubs) == args.top:
                _assign_tax()
                continue

            # stop when bit score drops below k fraction from top hit
            if args.top < 1 and (cbit - float(x[11])) / cbit > args.top:
                _assign_tax()
                continue

            # add current subject to collection
            csubs.append(x[1])

        else:
            # total hits below threshold (natural stop)
            if not skip:
                _assign_tax()

            # record top hit of new query
            cqry, csubs, cbit, skip = x[0], [x[1]], float(x[11]), False

    # handle last query
    if not skip:
        _assign_tax()


if __name__ == "__main__":
    main()
