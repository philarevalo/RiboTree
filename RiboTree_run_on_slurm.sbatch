#!/bin/bash

#SBATCH -p {name of your partition here}
#SBATCH -n 1

maxjobs=20
conda env create -f RiboTree_Environment.yml
source activate RiboTree
snakemake --snakefile RiboTree_IdentifyRibo.Snakefile --cluster-config cluster.yml --cluster "sbatch -p {cluster.partition} -N {cluster.nodes} -n {cluster.cores} --mem={cluster.mem}" --jobname {rulename}.{jobid} --jobs ${maxjobs}
snakemake --snakefile RiboTree_BuildPhylo.Snakefile --cluster-config cluster.yml --cluster "sbatch -p {cluster.partition} -N {cluster.nodes} -n {cluster.cores} --mem={cluster.mem}" --jobname {rulename}.{jobid} --jobs ${maxjobs}


