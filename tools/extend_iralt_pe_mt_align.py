#!/usr/bin/env python
'''
    Format of PE-MT alignments from IR-ALT are not compatible with qe-corpus-builder programs. In detail:
    -- Null-aligned tokens are not explicitly written in IR-ALT format.
    -- PE and MT indices are not arranged in ascending order respectively.

    In order to keep compatible, this script do the following things:
    --  Re-arrange the alignments, keeping MT indices in ascending order and adding Null-MT alignments.
    --  For PE-Null alignments, technically they can be inserted into the next position of (PE-1)-* or
        the previous position of (PE+1)-*. Currently, we adopt the first strategy.
'''

import argparse
import collections

def search(pairs, target_pe):
    '''
    Because order of PE indices is not assured in pairs.
    So we can only do a linear search rather than binary search.
    '''
    for i, (pe, mt) in enumerate(pairs):
        if pe == target_pe:
            return i
    return -1

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('-a', '--align',
                        help="Path to the original PE-MT alignment file.")
    parser.add_argument('-p', '--pe',
                        help="Path to the PE corpus file.")
    parser.add_argument('-m', '--mt',
                        help="Path to the MT corpus file.")
    parser.add_argument('-o', '--output',
                        help="Path to the output file (extended PE-MT alignments compatible with QCB).")

    args = parser.parse_args()
    return args

def main():
    def read_file(fn):
        with open(fn, 'r') as f:
            return [l.strip() for l in f]

    args = parse_args()

    pe_mt_lines = read_file(args.align)
    pe_lines = read_file(args.pe)
    mt_lines = read_file(args.mt)

    wf = open(args.output, 'w')

    for line_no in range(len(pe_mt_lines)):
        pe2mt_align_dict = collections.defaultdict(list)
        mt2pe_align_dict = collections.defaultdict(list)
        for align in pe_mt_lines[line_no].strip().split():
            a, b = map(int, align.split('-'))
            pe2mt_align_dict[a].append(b)
            mt2pe_align_dict[b].append(a)

        pe_token_len = len(pe_lines[line_no].strip().split())
        mt_token_len = len(mt_lines[line_no].strip().split())
        align_pairs = []
        for mt in range(mt_token_len):
            if len(mt2pe_align_dict[mt]) > 0:
                align_pairs.extend([(pe, mt) for pe in sorted(mt2pe_align_dict[mt])])
            else:
                align_pairs.append((None, mt))
        for pe in range(pe_token_len):
            if len(pe2mt_align_dict[pe]) == 0:
                if pe == 0: pos = 0
                else: pos = search(align_pairs, pe - 1) + 1
                align_pairs.insert(pos, (pe, None))

        align_pair_strs = []
        for a, b in align_pairs:
            sa = f"{a}" if a is not None else ''
            sb = f'{b}' if b is not None else ''
            align_pair_strs.append(f"{sa}-{sb}")
        wf.write(' '.join(align_pair_strs) + '\n')

    wf.close()

if __name__ == '__main__':
    main()