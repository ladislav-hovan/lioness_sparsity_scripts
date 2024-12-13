#!/usr/bin/env python

### Imports ###
import numpy as np
import pandas as pd

from collections import Counter
from pathlib import Path
from typing import Callable, Optional, Sequence

### Functions ###
def resample_reads_control(
    series: pd.Series,
    rng: np.random.Generator,
) -> pd.Series:
    """
    Generates a control expression pandas Series by sampling from the
    distribution of genes until the same number of reads is reached.
    Currently too slow to be used in practice.

    Parameters
    ----------
    series : pd.Series
        pandas Series with gene expressions
    rng : np.random.Generator
        numpy random Generator object to be used

    Returns
    -------
    pd.Series
        pandas Series with control gene expression
    """

    target_n = series.sum()
    probs = series.values / target_n

    selection = rng.choice(series.index, target_n, p=probs)
    counts = Counter(selection)

    res_series = pd.Series(counts.values(), index=counts.keys()
        ).reindex_like(series).fillna(0).astype(int)

    return res_series


def adjust_reads_control(
    series: pd.Series,
    rng: np.random.Generator,
    proportion: float = 0.1,
) -> pd.Series:
    """
    Generates a control expression pandas Series by perturbing the
    number of reads through multiplying by a number drawn from a normal
    distribution.

    Parameters
    ----------
    series : pd.Series
        pandas Series with gene expressions
    rng : np.random.Generator
        numpy random Generator object to be used
    proportion : float, optional
        Standard deviation for the normal distribution that will be used
        to perturb the gene expression data, by default 0.1

    Returns
    -------
    pd.Series
        pandas Series with control gene expression
    """

    offsets = series.apply(lambda x: round(rng.normal(0, x * proportion)))
    adj_series = series + offsets
    adj_series = adj_series.apply(lambda x: x if x > 0 else 1)

    return adj_series


def generate_control_expression(
    expression: Path,
    output: Sequence[Path],
    n_repeats: int,
    random_seed: Optional[int] = None,
    pert_function: Callable = adjust_reads_control,
) -> None:
    """
    Generates a control gene expression for multiple samples in multiple
    replicates by perturbing the given distribution of genes.

    Parameters
    ----------
    expression : Path
        Path to the file containing the expression data
    output : Sequence[Path]
        Paths to the files where the control expression data will be
        saved
    n_repeats : int
        Number of repeats to generate
    random_seed : Optional[int], optional
        Seed for the random number generator or None to use the system
        entropy, by default None
    pert_function : Callable, optional
        Function to be used to generate control gene expression,
        by default adjust_reads_control
    """

    rng = np.random.default_rng(random_seed)

    df = pd.read_csv(expression, sep='\t', index_col=0, header=None)

    for attempt in range(n_repeats):
        temp = {}
        for pos,col in enumerate(df.columns):
            temp[pos] = pert_function(df[col], rng)
        df_temp = pd.DataFrame({i: temp[i] for i in range(df.shape[1])})
        df_temp.to_csv(output[attempt], sep='\t', header=False)

### Main body ###
if __name__ == '__main__':
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument('-e', '--expression', dest='expression',
        help='file to load the gene expression from', metavar='FILE')
    parser.add_argument('-o', '--output', dest='output', nargs='+',
        help='files to save the control gene expression into', metavar='FILE')
    parser.add_argument('-nr', '--n-repeats', dest='n_repeats', type=int,
        help='number of repeats', metavar='INT')
    parser.add_argument('-rs', '--random-seed', dest='random_seed', type=int,
        help='random seed', metavar='INT', default=None)

    args = parser.parse_args()

    generate_control_expression(
        args.expression,
        args.output,
        args.n_repeats,
        args.random_seed,
    )