"""This robot is used to blastn the input CDS fasta file to the assembly database,
then get the location and other hsp information,
at last, extract the upstream CDS sequence information.--Xiao fei 06-04-2019"""

import sys
import argparse
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio import SeqRecord
from Bio import Seq
import copy
import subprocess


def parse_args():
    parser=argparse.ArgumentParser(description="Welcome to use Xiao_Fei_Robot")
    parser.add_argument('--CDS', required=True, type=str, metavar='FILENAME', help="the CDS fasta filename you want to blast")
    parser.add_argument('--DB', required=True, type=str, metavar='FILENAME', help="the database name in the command blastn -db")
    parser.add_argument('--WASH_LIST', required=True, type=str, metavar='FILENAME', help="the filename list you want to wash and get the upstream of CDS only")
    parser.add_argument('--CUT', default=100, type=int, metavar='up-stream cut length', help="the length(bp) you want to cut upstream of the CDS")
    parser.add_argument('--OUT', default="xiao_robot_extract_beside_CDS", type=str, metavar='directory', help="Output directory name")
    return parser.parse_args()


def main():
    args=parse_args()
    subprocess.run(["mkdir", "./"+args.OUT], check=True)
#doing blast
    do_blastn=NcbiblastnCommandline(query=args.CDS, outfmt=5, db=args.DB, out="temp_blastout.xml")
    stdout, stderr = do_blastn()
    if stderr:
        print("blastn meet some error")
    else:
        print("blastn step passed")

#extract information from blast.xml output lists of result [title, cutlocation1, cutlocation2, CDSname]    
    extracted_results=[]
    temp_blast_handle = open("temp_blastout.xml")
    blast_records = NCBIXML.parse(temp_blast_handle)
    #temp_blast_handle.close()
    for blast_record in blast_records:
        # Do something with blast_record
        for alignment in blast_record.alignments:
            title=alignment.title.split(' ',1)[1]
            for hsp in alignment.hsps:
                if hsp.identities / hsp.align_length >= 0.90:
                    if hsp.sbjct_start - hsp.sbjct_end <=0:
                        extracted_result=[title, hsp.sbjct_start-args.CUT, hsp.sbjct_start, blast_record.query, 0] #remeber this is 1-base system
                        extracted_results.append(extracted_result)
                    else:
                        extracted_result=[title, hsp.sbjct_start, hsp.sbjct_start+args.CUT, blast_record.query, 1] #remeber this is 1-base system
                        extracted_results.append(extracted_result)
                    break # this means if it get one hsp larger than 0.95 than out of the for alignment.hsps loop avoid two or more hsp with same id
                else:
                    pass
    print("extract information from blast.xml output lists of result [title, cutlocation1, cutlocation2, CDSname, strands(0 or 1)] passed")
    #print(extracted_results)
    
# using the extracted_results to extract the seq information
    wash_list_handle = open(args.WASH_LIST)
    wash_list = wash_list_handle.read().splitlines()
    #wash_list_handle.close()
    for wash_name in wash_list:
        new_seq_records=[]
        seq_records=SeqIO.parse(wash_name,"fasta")
        for seq_record in seq_records:
            for extracting in extracted_results:
                if extracting[0] == seq_record.description:
                    new_seq_record = copy.deepcopy(seq_record)
                    #new_seq_record.description = seq_record.description + "---FIX---" + extracting[3]
                    new_seq_record.description = wash_name + "---FIX---" + extracting[3]
                    new_seq_record.seq = seq_record.seq[int(extracting[1])-1 : int(extracting[2])-1] #remeber this is 0-base system
                    if extracting[4] == 1:
                        new_seq_record.seq = new_seq_record.seq.reverse_complement()
                    if len(new_seq_record.seq) == args.CUT:  #clean the ones which cannot fully cut with the cutting length, this can be happened at the two ends of one contig
                        new_seq_records.append(new_seq_record)
                    else:
                        pass
                else:
                    pass
        SeqIO.write(new_seq_records,"washed_"+wash_name, "fasta")
        subprocess.run(["mv", "washed_"+wash_name, "./"+args.OUT], check=True)
    print("dear xiao, you washing step is finished version 1.1, enjoy!")
   
if __name__ == '__main__':
    sys.exit(main())
   









