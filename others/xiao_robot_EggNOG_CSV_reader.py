'''this robot is used to deal with the data transformation of the output of EggNOG annotation file'''

import sys
import csv
import argparse


def parse_args():
    parser=argparse.ArgumentParser(description="Welcome to use Xiao_Fei_Robot")
    parser.add_argument('--input', required=True, type=str, metavar='FILENAME', help="the EggNOG output annotation file")
    parser.add_argument('--GET_KEGG_ID', action='store_const', const=True, metavar='SAVE_SOME_TEMP', help="get the KEGG ids")
    parser.add_argument('--GET_WEGO_ID', action='store_const', const=True, metavar='SAVE_SOME_TEMP', help="get the GO terms according name_id for the down stream of wego online")
    parser.add_argument('--GET_WEGO_CLEAN_ID', action='store_const', const=True, metavar='SAVE_SOME_TEMP', help="get the GO terms according name_id for the down stream of wego online, clean means filter the query without hit the GO term")
    #parser.add_argument('--WASH_LIST', required=True, type=str, metavar='FILENAME', help="the filename list you want to wash and get the upstream of CDS only")
    #parser.add_argument('--CUT', default=100, type=int, metavar='up-stream cut length', help="the length(bp) you want to cut upstream of the CDS")
    #parser.add_argument('--OUT', default="xiao_robot_extract_beside_CDS", type=str, metavar='directory', help="Output directory name")
    return parser.parse_args()


def get_kegg_id (inputfile_name):
    with open (inputfile_name, newline='') as csvfile:
        Fieldnames= [
            'query', 
            'Seed_Orthology', 
            'evalue', 
            'score', 
            'Predicted_name', 
            'GO_terms', 
            'KEGG_KO', 
            'BiGG_reactions', 
            'tax_scope', 
            'eggNOG_OGs', 
            'best_OG', 
            'COG_Cat', 
            'eggNOG_HMM_Desc']
        reader = csv.DictReader(csvfile, delimiter='\t', fieldnames= Fieldnames)
        for row in reader:
            KEGG_ID = row['KEGG_KO'].split(",")
            for id in KEGG_ID:
                print(id, file=open(inputfile_name+".KEGG_ID", "a") )

def get_query_link_GO_term (inputfile_name):
    with open (inputfile_name, newline='') as csvfile:
        Fieldnames= [
            'query', 
            'Seed_Orthology', 
            'evalue', 
            'score', 
            'Predicted_name', 
            'GO_terms', 
            'KEGG_KO', 
            'BiGG_reactions', 
            'tax_scope', 
            'eggNOG_OGs', 
            'best_OG', 
            'COG_Cat', 
            'eggNOG_HMM_Desc']
        reader = csv.DictReader(csvfile, delimiter='\t', fieldnames= Fieldnames)
        for row in reader:
            name_id = row['query']
            GO_TERMS = row['GO_terms'].split(",")
            for id in GO_TERMS:
                print(name_id + '\t' + id, file=open(inputfile_name+".wego", "a") )

def get_query_link_GO_term_clean (inputfile_name):
    with open (inputfile_name, newline='') as csvfile:
        Fieldnames= [
            'query', 
            'Seed_Orthology', 
            'evalue', 
            'score', 
            'Predicted_name', 
            'GO_terms', 
            'KEGG_KO', 
            'BiGG_reactions', 
            'tax_scope', 
            'eggNOG_OGs', 
            'best_OG', 
            'COG_Cat', 
            'eggNOG_HMM_Desc']
        reader = csv.DictReader(csvfile, delimiter='\t', fieldnames= Fieldnames)
        for row in reader:
            name_id = row['query']
            GO_TERMS = row['GO_terms'].split(",")
            #print(GO_TERMS)
            if GO_TERMS == ['']:
                pass
                #print("testinginging")
            else:
                for id in GO_TERMS:
                    print(name_id + '\t' + id, file=open(inputfile_name+".clean.wego", "a") )

            

def main():
    args=parse_args()
    if args.GET_KEGG_ID:
        get_kegg_id (args.input)
    if args.GET_WEGO_ID:
        get_query_link_GO_term (args.input)
    if args.GET_WEGO_CLEAN_ID:
        get_query_link_GO_term_clean (args.input)




print("hi xiao you job is finished! cheers!")
if __name__ == '__main__':
    sys.exit(main())

