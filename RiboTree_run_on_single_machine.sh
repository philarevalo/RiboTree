#!/bin/bash

conda env create -f RiboTree_Environment.yml
source activate RiboTree
snakemake --snakefile RiboTree_IdentifyRibo.Snakefile
snakemake --snakefile RiboTree_BuildPhylo.Snakefile


