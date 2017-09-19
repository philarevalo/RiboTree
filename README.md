RiboTree
=======================

This pipeline takes as input a directory of genomes in complete or draft format (one file = one genome) and outputs a phylogeny of those genomes based on single copy ribosomal proteins. Ribosomal proteins are identified using [hmmsearch](http://hmmer.org/) and an alignment of ribosomal proteins from [Yutin et al., 2012](http://dx.doi.org/10.1371/journal.pone.0036972), aligned using [mafft](http://mafft.cbrc.jp/alignment/software/), and used to build a phylogeny with [raxml](http://sco.h-its.org/exelixis/software.html). The basis for this pipeline is the method detailed in [Hehemann, et al., 2016](http://www.nature.com/articles/ncomms12860). It utilizes [snakemake](http://snakemake.readthedocs.io/en/latest/) for ease of use and reproducibility.

##How to run

Dependencies of this pipeline are handled by conda. You can download miniconda (for python 3.5 or higher please!) [here](https://conda.io/miniconda.html).

After installing miniconda:

1. Modify the `RiboTree_config.yml` file to reflect appropriate filepaths, parameters, and filenames. All parameters are described in comments within the `RiboTree_config.yml` file.

1b. If running on a cluster that uses slurm as a job scheduler, modify the `RiboTree_cluster_params.yml` file to match your cluster setup.

2a. If running on a single machine, simply run the `RiboTree_run_on_single_machine.sh` shell script.

2b. If running on a cluster that uses slurm as a job scheduler, modify the `RiboTree_run_on_slurm.sbatch` jobscript as appropriate and submit to the job scheduler.


##Notes

Not tested on Windows or OSX machines. I suspect that at present the code will break in windows due to the use of `\` as opposed to `/` in filepaths on Windows systems. May work in OSX. test.