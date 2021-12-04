#!/usr/bin/env python
import sys

def main():
    lines = sys.stdin.read().strip().split('\n')
    for l in lines:
        aligns = [tuple(map(int, a.split('-'))) for a in l.strip().split()]
        rev_aligns = [(a[1], a[0]) for a in aligns]
        rev_aligns.sort()
        print(' '.join([f'{a}-{b}' for a, b in rev_aligns]))

if __name__ == '__main__':
    main()