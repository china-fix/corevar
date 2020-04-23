'''this robot is modified from the script of xiao_robot_extract_beside_CDS.
in this situation, robot can help blastn and extract the CDS in each strain.
xiao fei 2019-04-14
''' 

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
    parser.add_argument('--WASH_LIST', required=True, type=str, metavar='FILENAME', help="the filename list you want to wash and get the CDS")
    #parser.add_argument('--CUT', default=100, type=int, metavar='up-stream cut length', help="the length(bp) you want to cut upstream of the CDS")
    parser.add_argument('--OUT', default="xiao_robot_extract_CDS", type=str, metavar='directory', help="Output directory name")
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

#extract information from blast.xml output lists of result    
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
                        extracted_result=[title, hsp.sbjct_start, hsp.sbjct_end, blast_record.query, 0] #remeber this is 1-base system
                        extracted_results.append(extracted_result)
                    else:
                        #print("tesingingingingingi")
                        extracted_result=[title, hsp.sbjct_end, hsp.sbjct_start, blast_record.query, 1] #remeber this is 1-base system
                        extracted_results.append(extracted_result)
                else:
                    pass
    print("extract information from blast.xml output lists of result [title, start_location, end_location, CDSname, strand(1 or 0)] passed")
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
                    new_seq_record.seq = seq_record.seq[int(extracting[1])-1 : int(extracting[2])] #remeber this is 0-base system
                    if extracting[4] == 1:
                        #print("testingingingingingi")
                        new_seq_record.seq = new_seq_record.seq.reverse_complement()
                    new_seq_records.append(new_seq_record)
                else:
                    pass
        SeqIO.write(new_seq_records,"washed_"+wash_name, "fasta")
        subprocess.run(["mv", "washed_"+wash_name, "./"+args.OUT], check=True)
    print("dear xiao, you washing step is finished version 1.2, enjoy!")
   
if __name__ == '__main__':
    sys.exit(main())
   







