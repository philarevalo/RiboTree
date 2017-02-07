from Bio import SeqIO
from collections import defaultdict
import argparse
import os

def main():
    parser = argparse.ArgumentParser(
        description=('Makes protein files from hmmer output'),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('-c', '--input_fasta', help='input fasta file of all ORFs')
    parser.add_argument('-m' , '--hmmer_file', help='hmmer output file')
    parser.add_argument('-o', '--out_path', help='output directory')
    parser.add_argument('-d', '--dummy_file', help='dummy file')
    args = parser.parse_args()

    infasta = args.input_fasta
    inhmm = args.hmmer_file
    outpath = args.out_path
    dummy_file = args.dummy_file
    seq_dict = {s.id : s for s in SeqIO.parse(infasta, 'fasta')}
    hmm_dict = defaultdict(list)


    with open(inhmm, 'r') as infile:
        for line in infile:
            if line[0] != '#':
                prot_name = line.split()[0]
                ribo_prot = line.split()[2]
                hmm_dict[ribo_prot].append(seq_dict[prot_name])
    all_strains = [s.id.split('_')[0] for prot, seqlist in hmm_dict.items() for s in seqlist]
    num_strains = len(set(all_strains))
    for prot, seqlist in hmm_dict.items():
        outname = outpath + prot + '_unaligned.fasta'
        final_seqlist = []
        for s in seqlist:
            strain = s.id.split('_')[0]
            s.id = strain
            final_seqlist.append(s)
        if len(final_seqlist) > 0.5 * num_strains:
            SeqIO.write(final_seqlist, outname, 'fasta')
    os.system('touch ' + dummy_file)

if __name__ == '__main__':
	main()