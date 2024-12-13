#!/usr/bin/env python

### Imports ###
import numpy as np
import scipy.stats as stats

from pathlib import Path
from typing import Iterable

from lib.functions import load_file

### Functions ###
def calculate_correlations(
    baseline: Path,
    to_compare: Iterable[Path],
    output_pearson: Path,
    output_spearman: Path,
) -> None:
    """
    Calculates the Pearson and Spearman correlations of the values in
    given files compared to the baseline.

    Parameters
    ----------
    baseline : Path
        Path to the file with baseline values
    to_compare : Iterable[Path]
        Iterable of Paths to the files to be compared
    output_pearson : Path
        Path to the file where Pearson correlations will be saved
    output_spearman : Path
        Path to the file where Spearman correlations will be saved
    """

    baseline_df = load_file(baseline)

    pearson_r = []
    spearman_r = []
    for file in to_compare:
        df = load_file(file).fillna(0)
        for col in df.columns:
            pearson_r.append(stats.pearsonr(baseline_df[col], df[col])[0])
            spearman_r.append(stats.spearmanr(baseline_df[col], df[col])[0])

    np.save(output_pearson, pearson_r)
    np.save(output_spearman, spearman_r)

### Main body ###
if __name__ == '__main__':
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument('-b', '--baseline', dest='baseline',
        help='file to compare the others to', metavar='FILE')
    parser.add_argument('-tc', '--to-compare', dest='to_compare', nargs='+',
        help='files to be compared to baseline', metavar='FILE')
    parser.add_argument('-op', '--output-pearson', dest='output_pearson',
        help='file to save the Pearson correlations to', metavar='FILE')
    parser.add_argument('-os', '--output-spearman', dest='output_spearman',
        help='file to save the Spearman correlations to', metavar='FILE')

    args = parser.parse_args()

    calculate_correlations(
        args.baseline,
        args.to_compare,
        args.output_pearson,
        args.output_spearman,
    )