#!/usr/bin/env python

### Imports ###
import pandas as pd

from pathlib import Path

### Functions ###
def mark_zero_expression_genes(
    expression: Path,
    output: Path,
) -> None:


    df = pd.read_csv(expression, sep='\t', index_col=0, header=None)

    is_zero = (df == 0)

    is_zero.to_csv(output, sep='\t', header=False)

### Main body ###
if __name__ == '__main__':
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument('-e', '--expression', dest='expression',
        help='file to load the gene expression from', metavar='FILE')
    parser.add_argument('-o', '--output', dest='output',
        help='file to save the marks into', metavar='FILE')

    args = parser.parse_args()

    mark_zero_expression_genes(
        args.expression,
        args.output,
    )