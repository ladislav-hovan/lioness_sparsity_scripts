#!/usr/bin/env python

### Imports ###
import numpy as np
import pandas as pd

from collections import Counter
from pathlib import Path
from typing import Optional, Sequence

### Functions ###
def resample_reads_single(
    series: pd.Series,
    sparsity: float,
    rng: np.random.Generator,
    min_chunk: int = 10000,
) -> pd.Series:
    """
    Sparsifies a pandas Series to a desired level by sampling from the
    distribution of genes until the sparsity level is reached.

    Parameters
    ----------
    series : pd.Series
        pandas Series with gene expression
    sparsity : float
        Desired level of sparsity in percentage points
    rng : np.random.Generator
        numpy random Generator object to be used
    min_chunk : int, optional
        minimal number of reads to generate in a single step,
        by default 10000

    Returns
    -------
    pd.Series
        pandas Series with sparsified gene expression
    """

    target_n = int((1.0 - sparsity) * len(series))
    probs = series.values / series.sum()
    counts = Counter()

    while len(counts) < target_n:
        to_gen = target_n - len(counts)
        selection = rng.choice(series.index, max(to_gen, min_chunk), p=probs)
        temp = counts.copy()
        temp.update(selection)
        if len(temp) < target_n:
            counts = temp
        else:
            index = 0
            while len(counts) < target_n:
                counts[selection[index]] += 1
                index += 1

    sparse_series = pd.Series(counts.values(), index=counts.keys()
        ).reindex_like(series).fillna(0).astype(int)

    return sparse_series


def resample_reads(
    expression: Path,
    output: Sequence[Path],
    sparsity: float,
    n_repeats: int,
    random_seed: Optional[int] = None,
) -> None:
    """
    Sparsifies gene expression for multiple samples in multiple
    replicates by sampling from the distribution of genes until the
    desired sparsity level is reached.

    Parameters
    ----------
    expression : Path
        Path to the file containing the expression data
    output : Sequence[Path]
        Paths to the files where the sparse expression data will be
        saved
    sparsity : float
        Desired level of sparsity in percentage points
    n_repeats : int
        Number of repeats to generate
    random_seed : Optional[int], optional
        Seed for the random number generator or None to use the system
        entropy, by default None
    """

    rng = np.random.default_rng(random_seed)

    df = pd.read_csv(expression, sep='\t', index_col=0, header=None)

    for attempt in range(n_repeats):
        temp = {}
        for pos,col in enumerate(df.columns):
            temp[pos] = resample_reads_single(df[col], 0.01 * sparsity, rng)
        df_temp = pd.DataFrame({i: temp[i] for i in range(df.shape[1])})
        df_temp.to_csv(output[attempt], sep='\t', header=False)

### Main body ###
if __name__ == '__main__':
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument('-e', '--expression', dest='expression',
        help='file to load the gene expression from', metavar='FILE')
    parser.add_argument('-o', '--output', dest='output', nargs='+',
        help='files to save the sparse gene expression into', metavar='FILE')
    parser.add_argument('-s', '--sparsity', dest='sparsity', type=float,
        help='level of sparsity in percentage points', metavar='FLOAT')
    parser.add_argument('-nr', '--n-repeats', dest='n_repeats', type=int,
        help='number of repeats', metavar='INT')
    parser.add_argument('-rs', '--random-seed', dest='random_seed', type=int,
        help='random seed', metavar='INT', default=None)

    args = parser.parse_args()

    resample_reads(
        args.expression,
        args.output,
        args.sparsity,
        args.n_repeats,
        args.random_seed,
    )