from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import glob
import pickle
from Bio.Alphabet import generic_dna
from collections import defaultdict
import argparse

def main():
    parser = argparse.ArgumentParser(
        description=('Concatenate a series of individually aligned genes into one alignment'),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('-i', '--input_prots', help='Space-delimited list of input aligned protein files')
    parser.add_argument('-m' , '--min_frac', help='Mininum fraction of column to be ungapped')
    parser.add_argument('-o', '--out_name', help='Output concatenated alignment filename')
    args = parser.parse_args()

    filelist = args.input_prots.split()
    min_frac = float(args.min_frac)
    out_name = args.out_name
    sequence_dict = defaultdict(str)
    final_seqs = []
    names1 = []
    all_strains = []
    num = 0

    for files in filelist:
        strains = [seqs.id for seqs in SeqIO.parse(files, 'fasta')]
        all_strains += strains
    all_strains = set(all_strains)    

    for files in filelist:
        prot_strains = []
        for seqs in SeqIO.parse(files, 'fasta'):
            taxon = seqs.id 
            prot_strains.append(taxon)
            sequence_dict[taxon] += str(seqs.seq)
            ali_length = len(seqs)
        gap_seq = "-" * ali_length
        for taxa in (all_strains - set(prot_strains)):
            sequence_dict[taxa] += gap_seq

    for name, sequence in sequence_dict.items():
        seq = Seq(sequence)
        print(len(seq))
        temp_record = SeqRecord(seq, id = name, description='')
        final_seqs.append(temp_record)

    SeqIO.write(final_seqs, out_name, 'fasta')

    columns_to_filter = []

    aln = AlignIO.read(out_name, 'fasta')
    for x in range(0, len(aln[0])):
        col = aln[ : , x]
        if col.count('-') >= min_frac * len(col):
            columns_to_filter.append(x)
    if len(columns_to_filter) > 0:
        new_align = aln[:, 0 : 0]

        for i, col in enumerate(sorted(columns_to_filter)):
            if i == 0:
                begin = 0
                new_align += aln[:, begin : col]
            elif i != len(columns_to_filter) - 1:
                begin = columns_to_filter[i - 1] + 1
                new_align += aln[:, begin : col]
            else:
                new_align += aln[:, col + 1 : ]
    else:
        new_align = aln
        
    AlignIO.write(new_align, out_name, 'fasta')
    
if __name__ == '__main__':
    main()