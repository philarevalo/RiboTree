"""
This snakemake pipeline takes as input a singleton hmmsearch result and outputs a
phylogeny (Raxml) built from an alignment of ribosomal proteins (mafft).
This technique was used to create a ribosomal reference tree 
in Hehemann et al., 2016 (doi:10.1038/ncomms12860).

"""

import glob
from collections import Counter

configfile:
    "RiboTree_config.yml"

all_prots = []
all_strains = []
with open(config["output_directory"] + "/hmmer_out/" + config["result_file_prefix"] + "_singletons.txt", 'r') as infile: 
    for line in infile:
        prot = line.split()[2]
        strain = line.split()[0].split('_')[0]
        all_strains.append(strain)
        all_prots.append(prot)
prot_counts = Counter(all_prots)
all_strains = set(all_strains)
all_prots = [prot for prot, count in prot_counts.items() if count > 0.5 * len(all_strains)]

rule target:    
    input:
        "RAxML_bestTree." + config['result_file_prefix'] + ".ribo.mafft.concat.tre"

rule run_raxml:
    input:
        config['result_file_prefix'] + ".ribo.mafft.concat.fasta"
    output:
        "RAxML_bestTree." + config['result_file_prefix'] + ".ribo.mafft.concat.tre"
    shell:
        "raxmlHPC {config[raxml_params]} -s {input} -n {output}"

rule concatenate_alignment:
    input:
        expand(config["output_directory"] + "/ribo_aligned/{prot_name}_aligned.fasta", prot_name = all_prots)
    output:
        config['result_file_prefix'] + ".ribo.mafft.concat.fasta"
    shell:
        "python {config[script_location]}/concatenate_alignment.py -m {config[min_col_frac]} -i '{input}' -o {output}"

rule align_prots:
    input:
        config["output_directory"] + "/ribo_unaligned/dummy_file.txt"
    output:
        config["output_directory"] + "/ribo_aligned/{prot_name}_aligned.fasta"
    run:
        p = str(output).split('/')[-1].split('_')[0]
        shell("mafft-linsi {config[output_directory]}/ribo_unaligned/{p}_unaligned.fasta > {output}")

rule make_ribo_files:
    input:
        singletons = config["output_directory"] + "/hmmer_out/" + config["result_file_prefix"] + "_singletons.txt",
        concatenated_orfs = config["output_directory"] + config["result_file_prefix"] + "_concatenated_orf.nt.fasta"
    output:
        config["output_directory"] + "/ribo_unaligned/dummy_file.txt"
    shell:
        "python {config[script_location]}/make_ribo_files.py -m {input[singletons]} -c {input[concatenated_orfs]} -o {config[output_directory]}/ribo_unaligned/ -d {output} " 
