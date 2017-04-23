from Bio import SeqIO
from collections import defaultdict, Counter
import numpy as np
import argparse
import os

def main():
    parser = argparse.ArgumentParser(
        description=('Finds singletons'),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('-i', '--input_hmm', help='input fasta file of all ORFs')
    parser.add_argument('-s', '--singleton_out_path', help='output file for singletons')
    args = parser.parse_args()
    cutoff = 1E-10
    factor = 100
    hmmer =  args.input_hmm
    temp_output = 'temp_output.txt' 
    check_for_multiple_hits(hmmer, factor, cutoff, temp_output)
    find_paralogs(temp_output, args.singleton_out_path)
    os.system('rm ' + temp_output)


def check_for_multiple_hits(in_hmm, factor_diff, min_eval, outfile_name):
    """
    Checks hmmsearch results for proteins that are assigned to multiple hmms.
    Assigns the protein to the top-scoring hmm if it passes a particular e-value 
    threshhold AND has an evalue that is a factor of <factor_diff> less than 
    the next best e-value.
    """
    hit_dict = defaultdict(list)
    rows = []
    rows_to_filter = []

    with open(in_hmm, 'r') as infile:
        for line in infile:
            if line[0] != '#':
                prot_name = line.split()[0]
                e_val = float(line.split()[4])
                if e_val <= min_eval:
                    hit_dict[prot_name].append(line)

    for prot, hits in hit_dict.items():
        if len(hits) > 1:
            rows_to_filter.append(hits)
        else:
            rows += hits

    """
    This section of the code does the actual filtering
    """       
    for row_group in rows_to_filter:
        sorted_rows = sorted(row_group, key = key_function) #sorts rows by e-value
        e_vals = [float(row.split()[4]) for row in sorted_rows] #list of all e-values
        if e_vals[1]/e_vals[0] >= factor_diff: #checks to see if best evalue is small enough compared to next best hit
            rows += [sorted_rows[0]]

    with open(outfile_name, 'w') as outfile:
        outfile.writelines(sorted(rows))

def find_paralogs(in_hmm, singleton_out):# , paralog_out):
    para_dict = defaultdict(list)
    single = []
    paralogs = []

    with open(in_hmm, 'r') as infile:
        for line in infile:
            if line[0] != '#':
                prot_name = line.split()[0]
                strain = '_'.join(prot_name.split('_')[0:-2])
                ribo_prot = line.split()[2]
                para_dict[(strain, ribo_prot)].append(line)

    for key, line_list in para_dict.items():
        if len(line_list) == 1:
            single += line_list
        else:
            paralogs += line_list
    with open(singleton_out, 'w') as outfile:
        outfile.writelines(sorted(single, key = key_function_ribo_prot))
        
def key_function(R):
    return float(R.split()[4])

def key_function_ribo_prot(R):
    return (R.split()[0].split('_')[0], R.split()[2])         

if __name__ == '__main__':
    main()