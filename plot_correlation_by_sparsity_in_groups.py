#!/usr/bin/env python

### Imports ###
import pandas as pd
import numpy as np

from itertools import pairwise
from pathlib import Path
from typing import Iterable, Literal

from lib.functions import plot_boxplots_in_groups

### Functions ###
def plot_correlation_by_sparsity_in_groups(
    input_files: Iterable[Path],
    sparsity_levels: Iterable[float],
    corr_name: Literal['pearson', 'spearman'],
    output: Path,
) -> None:
    """
    Plots the statistics for the given correlation type from the
    provided files as a boxplot.

    Parameters
    ----------
    input_files : Iterable[Path]
        Iterable of Paths to files containing the data
    sparsity_levels : Iterable[float]
        Iterable of sparsity levels corresponding to the files
    corr_name : Literal['pearson', 'spearman']
        Which type of correlation was used
    output : Path
        Path to the file where the plot will be saved
    """

    df = pd.DataFrame()
    if 0 == sparsity_levels[0]:
        df[0] = np.load(input_files[0])
        df[0] = df[0].fillna(df[0].mean())
    n = len(sparsity_levels) - 1
    for sparsity,(ft,ff) in zip(sparsity_levels[1:], zip(input_files[1:1+n], input_files[1+n:])):
        df[f'{sparsity}+'] = np.load(ft)
        df[f'{sparsity}+'] = df[f'{sparsity}+'].fillna(df[f'{sparsity}+'].mean())
        df[f'{sparsity}-'] = np.load(ff)
        df[f'{sparsity}-'] = df[f'{sparsity}-'].fillna(df[f'{sparsity}-'].mean())

    fig,_ = plot_boxplots_in_groups(df, corr_name=corr_name)
    fig.savefig(output, dpi=300, bbox_inches='tight')

### Main body ###
# if __name__ == '__main__':
#     from argparse import ArgumentParser

#     parser = ArgumentParser()
#     parser.add_argument('-i', '--input', dest='input_files', nargs='+',
#         help='files to plot the data from', metavar='FILE')
#     parser.add_argument('-sl', '--sparsity_levels', dest='s_levels', nargs='+',
#         help='levels of sparsity corresponding to the files', metavar='FLOAT',
#         type=float)
#     parser.add_argument('-cn', '--corr-name', dest='corr_name',
#         help='name of the correlation to use', metavar='NAME',
#         choices=['pearson', 'spearman'])
#     parser.add_argument('-o', '--output', dest='output',
#         help='file to save the plot to', metavar='FILE')

#     args = parser.parse_args()

#     plot_correlation_by_sparsity_in_groups(
#         args.input_files,
#         args.s_levels,
#         args.corr_name,
#         args.output,
#     )