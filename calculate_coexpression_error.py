#!/usr/bin/env python

### Imports ###
import numpy as np

from pathlib import Path
from typing import Iterable

from lib.functions import load_file

### Functions ###
def calculate_coexpression_error(
    baseline: Path,
    to_compare: Iterable[Path],
    output: Path,
) -> None:
    """
    Calculates the absolute error in gene coexpression values
    compared to the baseline starting from expression files.

    Parameters
    ----------
    baseline : Path
        Path to the expression file with baseline expression values
    to_compare : Iterable[Path]
        Iterable of Paths to the expression files to be compared
    output : Path
        Path to the file where coexpression error will be saved
    """

    baseline_df = load_file(baseline)
    bl_coexpr_df = baseline_df.T.corr()

    coexpr_err = []
    for file in to_compare:
        df = load_file(file)
        coexpr_df = df.T.corr()
        diff_df = coexpr_df - bl_coexpr_df
        for pos,row in enumerate(diff_df.columns):
            for col in diff_df.columns[pos+1:]:
                coexpr_err.append(abs(diff_df.loc[row,col]))

    np.save(output, coexpr_err)

### Main body ###
if __name__ == '__main__':
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument('-b', '--baseline', dest='baseline',
        help='file to compare the others to', metavar='FILE')
    parser.add_argument('-tc', '--to-compare', dest='to_compare', nargs='+',
        help='files to be compared to baseline', metavar='FILE')
    parser.add_argument('-o', '--output', dest='output',
        help='file to save the coexpression error to', metavar='FILE')

    args = parser.parse_args()

    calculate_coexpression_error(
        args.baseline,
        args.to_compare,
        args.output,
    )