#!/usr/bin/env python

### Imports ###
import numpy as np
import scipy.stats as stats

from pathlib import Path
from typing import Iterable

from lib.functions import load_file

### Functions ###
def calculate_correlations_in_groups(
    baseline: Path,
    to_compare: Iterable[Path],
    mark_files: Iterable[Path],
    output_pearson_true: Path,
    output_pearson_false: Path,
    output_spearman_true: Path,
    output_spearman_false: Path,
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

    pearson_rt = []
    pearson_rf = []
    spearman_rt = []
    spearman_rf = []
    for file,mark in zip(to_compare, mark_files):
        mark_df = load_file(mark)
        mark_df.columns = mark_df.columns.astype(str)
        df = load_file(file).fillna(0)
        for col in df.columns:
            bl_t = baseline_df[col][mark_df[col]]
            df_t = df[col][mark_df[col]]
            pearson_rt.append(stats.pearsonr(bl_t, df_t)[0])
            spearman_rt.append(stats.spearmanr(bl_t, df_t)[0])
            bl_f = baseline_df[col][~mark_df[col]]
            df_f = df[col][~mark_df[col]]
            pearson_rf.append(stats.pearsonr(bl_f, df_f)[0])
            spearman_rf.append(stats.spearmanr(bl_f, df_f)[0])

    np.save(output_pearson_true, pearson_rt)
    np.save(output_pearson_false, pearson_rf)
    np.save(output_spearman_true, spearman_rt)
    np.save(output_spearman_false, spearman_rf)

### Main body ###
if __name__ == '__main__':
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument('-b', '--baseline', dest='baseline',
        help='file to compare the others to', metavar='FILE')
    parser.add_argument('-tc', '--to-compare', dest='to_compare', nargs='+',
        help='files to be compared to baseline', metavar='FILE')
    parser.add_argument('-mf', '--mark-files', dest='mark_files', nargs='+',
        help='files containing the marks', metavar='FILE')
    parser.add_argument('-opt', '--output-pearson-true', dest='output_pearson_true',
        help='file to save the Pearson correlations for true marks to', metavar='FILE')
    parser.add_argument('-opf', '--output-pearson-false', dest='output_pearson_false',
        help='file to save the Pearson correlations for false marks to', metavar='FILE')
    parser.add_argument('-ost', '--output-spearman-true', dest='output_spearman_true',
        help='file to save the Spearman correlations for true marks to', metavar='FILE')
    parser.add_argument('-osf', '--output-spearman-false', dest='output_spearman_false',
        help='file to save the Spearman correlations for false marks to', metavar='FILE')

    args = parser.parse_args()

    calculate_correlations_in_groups(
        args.baseline,
        args.to_compare,
        args.mark_files,
        args.output_pearson_true,
        args.output_pearson_false,
        args.output_spearman_true,
        args.output_spearman_false,
    )