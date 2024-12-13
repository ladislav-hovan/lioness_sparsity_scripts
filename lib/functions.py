### Imports ###
import glob
import os
import time

import matplotlib.pyplot as plt
import pandas as pd

from netZooPy.lioness import Lioness
from netZooPy.panda import Panda
from pathlib import Path
from typing import Literal, Optional, Tuple

### Functions ###
def load_file(
    filename: Path,
) -> pd.DataFrame:
    """
    Reads a tsv or feather file and tries to process it based on the
    detected column names.

    Parameters
    ----------
    filename : Path
        Path to the file of interest

    Returns
    -------
    pd.DataFrame
        Processed pandas DataFrame

    Raises
    ------
    ValueError
        If the extension is neither 'tsv' nor 'feather'
    """

    extension = filename.split('.')[-1]

    if extension == 'tsv':
        # Expression files
        df = pd.read_csv(filename, sep='\t', index_col=0, header=None)
    elif extension == 'feather':
        raw_df = pd.read_feather(filename)
        if 'gene1' in raw_df.columns:
            # Correlation networks
            df = raw_df.set_index(['gene1', 'gene2'])
        elif 'tf' in raw_df.columns and 'gene' in raw_df.columns:
            # PANDA networks
            df = raw_df.set_index(['tf', 'gene'])
        elif 'gene' in raw_df.columns:
            # Indegree files
            df = raw_df.set_index('gene')
        elif 'tf' in raw_df.columns:
            # Outdegree files
            df = raw_df.set_index('tf')
        else:
            # Assumed to be coexpression matrices
            df = raw_df.set_index('0')
    else:
        raise ValueError(f'Invalid extension for: {filename}, '
            'expected tsv or feather')

    return df


def get_most_recent_log_time(
) -> float:
    """
    Parses the name of the most recently modified snakemake log file
    in order to retrieve its creation time.

    Returns
    -------
    float
        Time indicated in the most recent log file in epoch seconds
    """

    # Get the most recently modified log file
    last_log = max(glob.glob(os.path.join('.snakemake', 'log', '*.log')),
        key=os.path.getmtime)
    # Parse the name to get the datetime string
    last_time = last_log.split('/')[-1].split('.')[0]
    # Convert the string to epoch time
    last_time_f = time.mktime(time.strptime(last_time, '%Y-%m-%dT%H%M%S'))

    return last_time_f


def create_lioness_networks(
    expression: pd.DataFrame,
    motif_prior: pd.DataFrame,
    ppi_prior: pd.DataFrame,
    output_dir: Path,
    computing: Literal['cpu', 'gpu'],
    lioness_options: dict,
) -> None:
    """
    Calculates the LIONESS networks from the provided input files using
    the provided settings and saves them in the provided output
    directory.

    Parameters
    ----------
    expression : pd.DataFrame
        pandas DataFrame containing the expression data
    motif_prior : pd.DataFrame
        pandas DataFrame containing the gene regulation prior
    ppi_prior : pd.DataFrame
        pandas DataFrame containing the PPI prior
    output_dir : Path
        Path to a directory where the output file will be saved
    computing : Literal["cpu", "gpu"]
        Which type of computing will be used
    lioness_options : dict
        Keyword options passed to the Lioness object
    """

    panda_obj = Panda(expression, motif_prior, ppi_prior, computing=computing,
        save_memory=False, keep_expression_matrix=True)
    _ = Lioness(panda_obj, computing=computing, save_dir=output_dir,
        save_fmt='npy', **lioness_options)


def plot_boxplots(
    data: pd.DataFrame,
    corr_name: Optional[Literal['pearson', 'spearman']],
    ylim: Tuple[float, float] = (0, 1),
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Plots the boxplot from the given data.

    Parameters
    ----------
    data : pd.DataFrame
        pandas DataFrame containing the data
    corr_name : Optional[Literal['pearson', 'spearman']]
        Which type of correlation was used to generate the data or None
        it it is coexpression error
    ylim : Tuple[float, float], optional
        Limits of the y axis, by default (0, 1)

    Returns
    -------
    Tuple[plt.Figure, plt.Axes]
        matplotlib Figure and Axes of the generated plot
    """

    fig,ax = plt.subplots(figsize=(8,4))

    ax.boxplot(data, labels=data.columns)
    ax.set_ylim(ylim)
    ax.set_xlabel('Sparsity % of expressed genes')
    if corr_name is not None:
        corr_name = corr_name[0].upper() + corr_name[1:]
        ax.set_ylabel(f'{corr_name} R compared to no sparsity')
    else:
        ax.set_ylabel('Abs. coexpr. error compared to no sparsity')
    ax.grid(axis='y')

    return fig,ax


def plot_boxplots_in_groups(
    data: pd.DataFrame,
    corr_name: Optional[Literal['pearson', 'spearman']],
    ylim: Tuple[float, float] = (0, 1),
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Plots the boxplot from the given data.

    Parameters
    ----------
    data : pd.DataFrame
        pandas DataFrame containing the data
    corr_name : Optional[Literal['pearson', 'spearman']]
        Which type of correlation was used to generate the data or None
        it it is coexpression error
    ylim : Tuple[float, float], optional
        Limits of the y axis, by default (0, 1)

    Returns
    -------
    Tuple[plt.Figure, plt.Axes]
        matplotlib Figure and Axes of the generated plot
    """

    fig,ax = plt.subplots(figsize=(8,4))

    # ax.boxplot(data[[data.columns[0]]], labels=data.columns[0],
    #     manage_ticks=False, positions=[0])
    ax.boxplot(data, positions=[0] + sum([[i-0.15, i+0.15] for i in range(1, 1 + (len(data.columns) - 1) // 2)], start=[]),
        manage_ticks=False, widths=0.3)
    ax.set_xticks(range((len(data.columns) + 1) // 2))
    ax.set_xticklabels([0] + [int(float(x[:-1])) for x in data.columns[1::2]])
    ax.set_ylim(ylim)
    ax.set_xlabel('Sparsity % of expressed genes')
    if corr_name is not None:
        corr_name = corr_name[0].upper() + corr_name[1:]
        ax.set_ylabel(f'{corr_name} R compared to no sparsity')
    else:
        ax.set_ylabel('Abs. coexpr. error compared to no sparsity')
    ax.grid(axis='y')

    return fig,ax