"""
This snakemake pipeline takes as input genomes in fasta format and outputs a filtered
hmmsearch result that indicates all single-copy ribosomal proteins found. The hmmsearch
is conducted against an hmm constructed from an alignment of ribosomal proteins first
reported in Yutin et al., 2012 (http://dx.doi.org/10.1371/journal.pone.0036972). This
technique was used to extract ribosomal proteins for a reference tree in Hehemann et al., 
2016 (doi:10.1038/ncomms12860).

N.B.!!! THIS EXPECTS GENOME FILES TO BE NAMED SUCH THAT THE FIRST WORD IN THE FILENAME
AS SPLIT BY UNDERSCORES IS THE STRAIN NAME, e.g.:

FF50.fasta
FF50_contigs.fasta

It also expects that all genome filenames follow the same formatting

"""
import glob
from Bio import SeqIO
import os

configfile:
    "RiboTree_config.yml"

contig_wildcard_string = config["contig_directory"] + "/*" + config["contig_extension"]
os.system('cp {infile} {outfile}'.format(
    infile = config['outgroup_genome'], 
    outfile = config["contig_directory"] + '/outgroup_genome' + config['contig_extension']))
strains = ['.'.join(f.split('/')[-1].split('.')[0:-1]) for f in glob.glob(contig_wildcard_string)]


rule target:
    input:
        config["output_directory"] + "/hmmer_out/" + config["result_file_prefix"] + "_singletons.txt"

rule filter_hmm_results:
    input:
        config["output_directory"] + "/hmmer_out/" + config["result_file_prefix"] + "_hmmer_result.tbl"
    output:
        config["output_directory"] + "/hmmer_out/" + config["result_file_prefix"] + "_singletons.txt"
    shell:
        "python {config[script_location]}/check_hmmsearch.py -i {input} -s {output}"

rule hmmsearch:
    input:
        concatenated_orfs = config["output_directory"] + config["result_file_prefix"] + "_concatenated_orf.aa.fasta",
        hmmfile = {config["ribo_hmm_file"]}
    output:
        text = config["output_directory"] + "/hmmer_out/" + config["result_file_prefix"] + "_hmmer_result.txt",
        table = config["output_directory"] + "/hmmer_out/" + config["result_file_prefix"] + "_hmmer_result.tbl",
        domains = config["output_directory"] + "/hmmer_out/" + config["result_file_prefix"] + "_hmmer_result.dom.tbl"
    shell:
        "hmmsearch -o {output[text]} --tblout {output[table]} --domtblout {output[domains]} {input[hmmfile]} {input[concatenated_orfs]}"

rule make_concatenated_orf_file:
    input:
        aa = expand(config["output_directory"] + "/orfs/{strain}_aa.fasta", strain = strains),
        nt = expand(config["output_directory"] + "/orfs/{strain}_nt.fasta", strain = strains)
    output:
        aa = config["output_directory"] + config["result_file_prefix"] + "_concatenated_orf.aa.fasta",
        nt = config["output_directory"] + config["result_file_prefix"] + "_concatenated_orf.nt.fasta"
    shell: 
        "cat {input[aa]} > {output[aa]}; cat {input[nt]} > {output[nt]}"

rule run_prodigal:
    input:
        config["contig_directory"] + "/{strain}" + config["contig_extension"] + ".reformatted"
    output:
        aa = config["output_directory"] + "/orfs/{strain}_aa.fasta",
        nt = config["output_directory"] + "/orfs/{strain}_nt.fasta"
    shell:
        "prodigal -i {input} -a {output[aa]} -d {output[nt]} -q"

rule reformat_contigs:
    #This chunk of code insures that sequence names in genome files
    #follow the correct format.
    input:
        genome_file = config["contig_directory"] + "/{strain}" + config['contig_extension'],
    output:
        output_file = config["contig_directory"] + "/{strain}" + config["contig_extension"] + ".reformatted"
    run:
        new_s = []
        strain = '.'.join(input["genome_file"].split('/')[-1].split('.')[0:-1])
        for i, s in enumerate(SeqIO.parse(input["genome_file"], 'fasta')):
            if len(s) > 0:
                s.id = strain + '_' + str(i)
                new_s.append(s)
        SeqIO.write(new_s, output["output_file"], 'fasta')