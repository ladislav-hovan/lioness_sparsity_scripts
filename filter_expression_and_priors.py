#!/usr/bin/env python

### Imports ###
import pandas as pd

from pathlib import Path

from lib.functions import load_file

### Functions ###
def filter_expression_and_priors(
    expression_path: Path,
    motif_prior_path: Path,
    ppi_prior_path: Path,
    expression_output: Path,
    motif_prior_output: Path,
    ppi_prior_output: Path,
) -> None:
    """
    Filters the gene expression file to only include genes with non-zero
    reads everywhere and also ensures matching between the gene
    expression and the motif and PPI priors.

    Parameters
    ----------
    expression : Path
        Path to the file containing the expression data
    motif_prior : Path
        Path to the file containing the gene regulation prior
    ppi_prior : Path
        Path to the file containing the PPI prior
    expression_output : Path
        Path to the file where the expression data will be saved
    motif_prior_output : Path
        Path to the file where the gene regulation prior will be saved
    ppi_prior_output : Path
        Path to the file where the PPI prior will be saved
    """

    # Remove genes with zero expression
    df = load_file(expression_path)
    mask = ((df > 0).sum(axis=1) == df.shape[1])
    df_filt = df.loc[mask]

    # Ensure correspondence between motif prior and expression
    motif_prior = pd.read_csv(motif_prior_path, sep='\t',
        names=['tf', 'gene', 'edge'])
    df_genes = set(df_filt.index)
    motif_prior_genes = set(motif_prior['gene'].unique())
    common_genes = df_genes.intersection(motif_prior_genes)
    df_common = df_filt.loc[list(common_genes)]
    mask = motif_prior['gene'].isin(common_genes)
    motif_prior_common = motif_prior.loc[mask]

    # Ensure all TFs from PPI prior in the motif prior
    ppi_prior = pd.read_csv(ppi_prior_path, sep='\t',
        names=['tf1', 'tf2', 'edge'])
    all_tfs = set(ppi_prior['tf1']).union(set(ppi_prior['tf2']))
    motif_prior_tfs = set(motif_prior_common['tf'].unique())
    extra_tfs = all_tfs - motif_prior_tfs
    mask = (~ppi_prior['tf1'].isin(extra_tfs) &
        ~ppi_prior['tf2'].isin(extra_tfs))
    ppi_prior_common = ppi_prior.loc[mask]

    # Save everything properly
    df_common.sort_index().to_csv(expression_output, sep='\t', header=False)
    motif_prior_common.sort_values(['tf', 'gene']).to_csv(motif_prior_output,
        sep='\t', index=False, header=False)
    ppi_prior_common.sort_values(['tf1', 'tf2']).to_csv(ppi_prior_output,
        sep='\t', index=False, header=False)

### Main body ###
if __name__ == '__main__':
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument('-e', '--expression', dest='expression',
        help='file to load the gene expression from', metavar='FILE')
    parser.add_argument('-mp', '--motif-prior', dest='motif_p',
        help='file to load the motif prior from', metavar='FILE')
    parser.add_argument('-pp', '--ppi-prior', dest='ppi_p',
        help='file to load the PPI prior from', metavar='FILE')
    parser.add_argument('-eo', '--expression-output', dest='expression_o',
        help='file to save the filtered gene expression into', metavar='FILE')
    parser.add_argument('-mpo', '--motif-prior-output', dest='motif_p_o',
        help='file to save the filtered motif prior into', metavar='FILE')
    parser.add_argument('-ppo', '--ppi-prior-output', dest='ppi_p_o',
        help='file to save the filtered PPI prior into', metavar='FILE')

    args = parser.parse_args()

    filter_expression_and_priors(
        args.expression,
        args.motif_p,
        args.ppi_p,
        args.expression_o,
        args.motif_p_o,
        args.ppi_p_o,
    )