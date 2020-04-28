import sys
import argparse
from Bio import Seq
from Bio import SeqIO
from Bio import SeqRecord
from Bio.Alphabet import generic_dna

"""This robot is used to extract user specified CDSs (a list file) from the Roary output of pangenome fasta file 
Xiao 29-03-2019"""

def parse_args():
    parser=argparse.ArgumentParser(description="Welcome to use Xiao_Fei_Robot")
    parser.add_argument('--pangenome', required=True, type=str, metavar='FILENAME', help="roary output pangenome fasta file")
    parser.add_argument('--extract', required=True, type=str, metavar='FILENAME', help="gene names you want to extract, each line one name")
    parser.add_argument('--output', required=True, type=str, metavar='FILENAME', help="Output filename")
    return parser.parse_args()



# triming the original pangenome fasta file
def main():
    args=parse_args()
    seq_records=SeqIO.parse(args.pangenome,"fasta")
    new_seq_records=[]
    for seq_record in seq_records:
        seq_record.id=seq_record.description.split()[1]
        new_seq_records.append(seq_record)
    SeqIO.write(new_seq_records,"temp_xiao_fei_robot", "fasta")

# get the dictionary of trimmed pangenome fasta
    seq_records_dict=SeqIO.index("temp_xiao_fei_robot","fasta")
    print(len(seq_records_dict))
    #print(list(seq_records_dict.keys()))

# extract the seq_records according the key_list
    #key_list=["sodC_1","group_3250"]
    with open(args.extract,'r') as f:
        key_list=[line.rstrip('\n') for line in f]
    extract_seq_records=[]
    extract_seq_records_translate=[]
    for key in key_list:
        extract_seq_records.append(seq_records_dict[key])
        temp=seq_records_dict[key]
        temp.seq=temp.seq.translate(table="Bacterial", to_stop=True)
        extract_seq_records_translate.append(temp)
    SeqIO.write(extract_seq_records,args.output,"fasta")
    SeqIO.write(extract_seq_records_translate,args.output+".p", "fasta")    


if __name__ == '__main__':
    sys.exit(main())
