#!/usr/bin/env python

### Imports ###
import pandas as pd

from pathlib import Path
from scipy.stats import halfnorm, uniform

### Functions ###
def add_expression_noise(
    expression: Path,
    output: Path,
    scale: float = 1.0,
) -> None:
    """
    Adds a normally distributed noise to the expression data.

    Parameters
    ----------
    expression : Path
        Path to the file containing the expression data
    output : Path
        Path to the file where the rescaled expression data will be
        saved
    sigma : float, optional
        Scale factor for the half-normal distribution of noise,
        by default 1.0
    """

    df = pd.read_csv(expression, sep='\t', index_col=0, header=None)

    # df_noise = halfnorm.rvs(scale=scale, size=df.shape)
    df_noise = uniform.rvs(scale=scale, size=df.shape)

    df_mod = df + df_noise

    df_mod.to_csv(output, sep='\t', header=False)

### Main body ###
if __name__ == '__main__':
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument('-e', '--expression', dest='expression',
        help='file to load the gene expression from', metavar='FILE')
    parser.add_argument('-o', '--output', dest='output',
        help='file to save the rescaled gene expression into', metavar='FILE')
    parser.add_argument('-s', '--scale', dest='scale', type=float,
        help='scale factor for the half-normal distribution',
        metavar='FLOAT', default=1.0)

    args = parser.parse_args()

    add_expression_noise(
        args.expression,
        args.output,
        args.scale,
    )