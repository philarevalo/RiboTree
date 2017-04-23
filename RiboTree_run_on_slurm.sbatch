#!/bin/bash

#SBATCH -p [insert partition name here]
#SBATCH -n 1

cluster_file=RiboTree_cluster_params.yml
maxjobs=20
conda env create -f RiboTree_Environment.yml
source activate ribo_env
snakemake --snakefile RiboTree_IdentifyRibo.Snakefile --cluster-config ${cluster_file} --cluster "sbatch -p {cluster.partition} -N {cluster.nodes} -n {cluster.cores} --mem={cluster.mem}" --jobname {rulename}.{jobid} --jobs ${maxjobs}
snakemake --snakefile RiboTree_BuildPhylo.Snakefile --cluster-config ${cluster_file} --cluster "sbatch -p {cluster.partition} -N {cluster.nodes} -n {cluster.cores} --mem={cluster.mem}" --jobname {rulename}.{jobid} --jobs ${maxjobs}