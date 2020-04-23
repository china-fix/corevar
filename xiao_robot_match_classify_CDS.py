'''
This robot can help to locally blastn the CDS.fasta to the reference genome and filter the orginal CDS.fasta to matched_CDS and un_matched_CDS  
'''

import subprocess
import sys
import argparse
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio import SeqRecord
from Bio import Seq
import copy


def parse_args():
    parser=argparse.ArgumentParser(description="Welcome to use Xiao_Fei_Robot")
    parser.add_argument('--CDS', required=True, type=str, metavar='FILENAME', help="the CDS fasta filename you want to blast")
    parser.add_argument('--REF', required=True, type=str, metavar='FILENAME', help="the complete reference genome of fasta format you want to blast to")
    parser.add_argument('--CUTOFF', default=0.9, type=float, metavar='DEFAULT 0.9', help="the lowest similarity value which classify as matched")
    parser.add_argument('--OUT', default="xiao_robot_match_classify_CDS", type=str, metavar='directory', help="Output directory name")
    return parser.parse_args()

#1. blast and get the xml
#2. parse the xml and get the matched and un_matched list
#3. extract and output the matched.fasta and un_matched.fasta

def doing_blast(query, reference):
    subprocess.run(["blastn", "-query", query, "-subject", reference, "-outfmt", "5", "-out", "temp_blast.xml"], check=True)


def filter_matching(cutoff):
    matched_list=[]
    temp_blast_handle = open("temp_blast.xml")
    blast_records = NCBIXML.parse(temp_blast_handle)
    for blast_record in blast_records:
        # Do something with blast_record
        for alignment in blast_record.alignments:
            hsp = alignment.hsps[0]
            if  hsp.identities / hsp.align_length >= cutoff and hsp.expect <= 1e-10 : #and hsp.query_start == 1 and hsp.query_end == blast_record.query_letters:
                matched_result=blast_record.query
                matched_list.append(matched_result)
            else:
                pass
    return matched_list



def main():
    args = parse_args()
    doing_blast(args.CDS, args.REF)
    print("##########################################")
    print("blasn and get the xml passed")

    matched_list = filter_matching(args.CUTOFF)
    print("##########################################")
    print("parse the xml and get the matched and un_matched list passed")

    
    seq_records=SeqIO.parse(args.CDS,"fasta")
    matched_seq_records=[]
    un_matched_seq_records=[]
    un_matched_list = []
    main_list=[]
    for seq_record in seq_records:
        main_list.append(seq_record.description)
        for matched_name in matched_list:
            if matched_name == seq_record.description:
                new_seq_record = copy.deepcopy(seq_record)
                matched_seq_records.append(new_seq_record)
    
    for main_name in main_list:
        if main_name not in matched_list:
            un_matched_list.append(main_name)     
    #print(un_matched_list)
    #print(seq_records)
    seq_records=SeqIO.parse(args.CDS,"fasta")
    for seq_record in seq_records:
        for un_matched_name in un_matched_list:
            if  un_matched_name == seq_record.description:
                new_seq_record = copy.deepcopy(seq_record)
                un_matched_seq_records.append(new_seq_record)

    
    SeqIO.write(matched_seq_records,"matched_" + args.OUT, "fasta")
    SeqIO.write(un_matched_seq_records,"un_matched_" + args.OUT, "fasta")
    print("##########################################")
    print("extract and output the matched and un_matched fasta file passed")

    print("hi xiao, you filtering work succeed! cheers!")


if __name__ == '__main__':
    sys.exit(main())
