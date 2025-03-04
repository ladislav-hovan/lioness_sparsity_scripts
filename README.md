# LIONESS sparsity testing scripts
These scripts are part of pipelines that measure the sensitivity of
LIONESS to varying levels of sparsity.


## Table of Contents
- [LIONESS sparsity testing scripts](#lioness-sparsity-testing-scripts)
  - [Table of Contents](#table-of-contents)
  - [General Information](#general-information)
  - [Features](#features)
  - [Setup](#setup)
  - [Usage](#usage)
  - [Project Status](#project-status)
  - [Room for Improvement](#room-for-improvement)
  - [Acknowledgements](#acknowledgements)
  - [Contact](#contact)
  - [License](#license)


## General Information
There are currently two pipelines using these scripts: one written in
[Snakemake](https://github.com/ladislav-hovan/lioness_sparsity_snakemake)
and one under construction in
[Nextflow](https://github.com/ladislav-hovan/lioness_sparsity_nextflow).
The pipelines test and summarise how the output of LIONESS changes 
when the gene expression data is modified to varying levels of sparsity. 
They are intended to automate this process fully, starting from input 
expression file, prior motif network, prior PPI (protein-protein 
interaction) network, and the settings provided.


## Features
The scripts perform the following tasks:
- Cleanup of expression data and prior networks
- Generation of sparsified expression
- Generation of LIONESS networks
- Summary Pearson and Spearman correlations comparing to the unmodified
  baseline:
  - edge weights
  - indegrees
  - expression
  - coexpression
- Average error in coexpression across sparsity levels


## Setup
The requirements are provided in a `requirements.txt` file.


## Usage
Please refer to the usage of the corresponding pipelines.
The scripts can also be used separately from these pipelines as they
implement a command line interface.
For more information on how to run them and the required input:

``` bash
# This particular script is only used as an example
python filter_expression_and_priors.py --help
```


## Project Status
The project is: _in progress_.


## Room for Improvement
Room for improvement:
- Finish refactoring shared functionality

To do:
- Expand the stratification of results to be more generally applicable


## Acknowledgements
Many thanks to the members of the 
[Kuijjer group](https://www.kuijjerlab.org/) 
at NCMM for their feedback and support.

This README is based on a template made by 
[@flynerdpl](https://www.flynerd.pl/).


## Contact
Created by Ladislav Hovan (ladislav.hovan@ncmm.uio.no).
Feel free to contact me!


## License
This project is open source and available under the 
[GNU General Public License v3](LICENSE).