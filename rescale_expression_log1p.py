#!/usr/bin/env python

### Imports ###
import numpy as np
import pandas as pd

from pathlib import Path

### Functions ###
def rescale_expression_log1p(
    expression: Path,
    output: Path,
) -> None:
    """
    Applies the log1p transformation to the provided expression data,
    which consists of adding 1 and then applying the natural logarithm.

    Parameters
    ----------
    expression : Path
        Path to the file containing the expression data
    output : Path
        Path to the file where the rescaled expression data will be
        saved
    """

    df = pd.read_csv(expression, sep='\t', index_col=0, header=None)

    df_mod = df.add(1).apply(np.log)

    df_mod.to_csv(output, sep='\t', header=False)

### Main body ###
if __name__ == '__main__':
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument('-e', '--expression', dest='expression',
        help='file to load the gene expression from', metavar='FILE')
    parser.add_argument('-o', '--output', dest='output',
        help='file to save the rescaled gene expression into', metavar='FILE')

    args = parser.parse_args()

    rescale_expression_log1p(
        args.expression,
        args.output,
    )