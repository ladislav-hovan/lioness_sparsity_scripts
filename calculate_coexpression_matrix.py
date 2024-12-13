#!/usr/bin/env python

### Imports ###
from pathlib import Path

from lib.functions import load_file

### Functions ###
def calculate_coexpression_matrix(
    expression: Path,
    output: Path,
) -> None:
    """
    Calculates the coexpression matrix from the given gene expression.

    Parameters
    ----------
    expression : Path
        Path to a file containing the gene expression
    output : Path
        Path where the coexpression matrix will be saved
    """

    df = load_file(expression)
    coexpr_df = df.T.corr()
    coexpr_df.reset_index().to_feather(output)

### Main body ###
if __name__ == '__main__':
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument('-e', '--expression', dest='expression',
        help='file to load the gene expression from', metavar='FILE')
    parser.add_argument('-o', '--output', dest='output',
        help='file to save the coexpression matrix into', metavar='FEATHER')

    args = parser.parse_args()

    calculate_coexpression_matrix(
        args.expression,
        args.output,
    )