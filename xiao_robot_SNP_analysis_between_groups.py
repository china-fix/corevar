"""this is a robot used to analysis the SNP between different strain groups. Tying to find the conserved SNP cases between the groups. 
1. extract specific CDS from strains and combine
2. Alignment with PRANK
3. SNP-sites
4. get score"""

import sys
import subprocess
import argparse
from Bio import SeqIO
import copy
import vcf
import os

def parse_args():
    parser=argparse.ArgumentParser(description="Welcome to use Xiao_Fei_Robot")
    parser.add_argument('--CDS_LIST', required=True, type=str, metavar='FILENAME', help="the CDS list you want to extract")
    parser.add_argument('--STRAIN_LIST', required=True, type=str, metavar='FILENAME', help="the strain list you want to extract from the fasta files (each file comtain all the CDS of a strain)")
    parser.add_argument('--GROUP_1', required=True, type=str, metavar='FILENAME', help="the strain (fasta files name) list you want to caculate as one group")
    parser.add_argument('--GROUP_2', required=True, type=str, metavar='FILENAME', help="the strain (fasta files name) list you want to caculate as one group")
    parser.add_argument('--TEMP_SAVE', action='store_const', const=True, metavar='SAVE_SOME_TEMP', help="this command help to save some temp infomation or files")
    parser.add_argument('--FAST', action='store_const', const=True, metavar='FAST ALIGNMENT WITH MAFFT', help="this command help to use mafft to do a fast alignment")
    #parser.add_argument('--CUT', default=100, type=int, metavar='up-stream cut length', help="the length(bp) you want to cut upstream of the CDS")
    parser.add_argument('--OUT', default="xiao_robot_SNP_analysis_between_groups", type=str, metavar='directory', help="Output directory name")
    return parser.parse_args()

def main():
    args=parse_args()
    subprocess.run(["mkdir", "./"+args.OUT], check=True)
    if args.TEMP_SAVE:
        subprocess.run(["mkdir", "./TEMP_SAVE"], check=True)

    #_1. extract specific CDS from strains and combine according different CDS
    with open(args.CDS_LIST) as cds_list:
        for cds in cds_list.read().splitlines():
            new_seq_records=[]
            with open(args.STRAIN_LIST) as strain_list:
                for strain in strain_list.read().splitlines():
                    seq_records=SeqIO.parse(strain,"fasta")
                    for seq_record in seq_records:
                        if cds == seq_record.description.split('---FIX---')[1]:
                            new_seq_record = copy.deepcopy(seq_record)
                            new_seq_record.id = seq_record.description.split(" ")[1]  # modifiy the fasta title for the downstream use
                            new_seq_record.description = ""                          #
                            new_seq_records.append(new_seq_record)
                        else:
                            pass
            SeqIO.write(new_seq_records,"temp_" + cds.split(' ')[0], "fasta")
            if args.TEMP_SAVE:
                subprocess.run(["cp", "temp_" + cds.split(' ')[0], "./TEMP_SAVE"], check=True)
    print("step 1. extract specific CDS from strains and combine according different CDS passed")

    #_2.Alignment with PRANK
    with open(args.CDS_LIST) as cds_list:
        for cds in cds_list.read().splitlines():
            if args.FAST:
                try:
                    #subprocess.run(["mafft",  "temp_" + cds.split(' ')[0] + ">" + "temp_" + cds.split(' ')[0]+".best.fas"], check=True) # 3.7 python can change to capture_output=True
                    os.system("mafft " + "temp_" + cds.split(' ')[0] + " > " + "temp_" + cds.split(' ')[0]+".best.fas")
                    if args.TEMP_SAVE:
                        subprocess.run(["cp", "temp_" + cds.split(' ')[0] + ".best.fas", "./TEMP_SAVE"], check=True)
                except subprocess.CalledProcessError:
                    pass
            else:
                try:
                    subprocess.run(["prank", "-d=" + "temp_" + cds.split(' ')[0], "-o=" + "temp_" + cds.split(' ')[0]], check=True) # 3.7 python can change to capture_output=True
                    if args.TEMP_SAVE:
                        subprocess.run(["cp", "temp_" + cds.split(' ')[0] + ".best.fas", "./TEMP_SAVE"], check=True)
                except subprocess.CalledProcessError:
                    pass
            subprocess.run(["rm", "temp_" + cds.split(' ')[0]], check=True)
            #subprocess.run(["mv", "temp_" + cds.split(' ')[0] + ".best.fas", "./"+args.OUT], check=True)
    print("#############################################################")
    print("step 2.Alignment with PRANK is passed")

    #_3.SNP calling with snp-sites, output the vcf files
    with open(args.CDS_LIST) as cds_list:
        vcf_filter_list=[]
        vcf_unfilter_list=[]
        vcf_unfilter_error_list=[]
        for cds in cds_list.read().splitlines():
            return_code = subprocess.run(["snp-sites", "-v" , "-otemp_" + cds.split(' ')[0] + ".vcf", "temp_" + cds.split(' ')[0] + ".best.fas"]).returncode
            if float(return_code) == 0:
                subprocess.run(["rm", "temp_" + cds.split(' ')[0] + ".best.fas"], check=True)
                if args.TEMP_SAVE:
                    subprocess.run(["cp", "temp_" + cds.split(' ')[0] + ".vcf", "./TEMP_SAVE"], check=True)
                subprocess.run(["mv", "temp_" + cds.split(' ')[0] + ".vcf", "./"+args.OUT], check=True)
                vcf_filter_list.append("temp_" + cds.split(' ')[0] + ".vcf")
            elif float(return_code) == 1:
                print("in "+ cds.split(' ')[0])
                try:
                    subprocess.run(["rm", "temp_" + cds.split(' ')[0] + ".best.fas"], check=True)
                    vcf_unfilter_list.append("temp_" + cds.split(' ')[0] + ".vcf")
                except subprocess.CalledProcessError:
                    vcf_unfilter_error_list.append("temp_" + cds.split(' ')[0] + ".vcf")
                    pass
                
            else:
                print(return_code)
                raise Exception("xiao robot meet some erros at the snp-sites command in step 3")
    '''for vcf_filter in vcf_filter_list:
        print(vcf_filter[5:-4], file=open("temp_vcf_filter_list", "a") )  
    for vcf_unfilter in vcf_unfilter_list:
        print(vcf_unfilter[5:-4], file=open("temp_vcf_unfilter_list", "a") ) '''
    print("#############################################################")
    print("step 3.SNP calling with snp-sites, output the vcf files passed")
    if args.TEMP_SAVE:
        with open("temp_vcf_filter_list", "w") as f:
            for s in vcf_filter_list:
                f.write(str(s) +"\n")
        with open("temp_vcf_unfilter_list", "w") as f:
            for s in vcf_unfilter_list:
                f.write(str(s) +"\n")
        with open("temp_vcf_unfilter_error_list", "w") as f:
            for s in vcf_unfilter_error_list:
                f.write(str(s) +"\n")
        subprocess.run(["mkdir", "./TEMP_SAVE/FILTER_LOG"], check=True)
        subprocess.run(["mv", "temp_vcf_filter_list", "./TEMP_SAVE/FILTER_LOG"], check=True)
        subprocess.run(["mv", "temp_vcf_unfilter_list", "./TEMP_SAVE/FILTER_LOG"], check=True)
        subprocess.run(["mv", "temp_vcf_unfilter_error_list", "./TEMP_SAVE/FILTER_LOG"], check=True)



    #_4.read and analysis the vcf files and get the scores
    # the function of this part is to get compare each SNP between each group and if the average difference is more than 0.75, add score as 1,then sum the scores in each CDS
    for vcf_filter in vcf_filter_list:
        with open("./"+args.OUT + "/"+ vcf_filter) as VCF_FILTER:
            vcf_reader = vcf.Reader(VCF_FILTER)
            analysis_score = 0
            zero_record = "###passed"
            for vcf_record in vcf_reader: 
                with open(args.GROUP_1) as group_1:
                    group_score_1 = 0
                    group_1_num = 0
                    error_1 = None
                    for name_1 in group_1.read().splitlines():
                        name_1 = name_1.split("ashed_")[1] 
                        name_1_follow = vcf_filter[5:-4]
                        try:
                            score_1 = vcf_record.genotype(name_1 + "---FIX---" + name_1_follow)["GT"]
                            group_score_1 += float(score_1)
                            group_1_num += 1
                        except KeyError as ve_1:
                            error_1 = ve_1
                            pass
                with open(args.GROUP_2) as group_2:
                    group_score_2 = 0
                    group_2_num = 0
                    error_2 = None
                    for name_2 in group_2.read().splitlines():
                        name_2 = name_2.split("ashed_")[1] 
                        name_2_follow = vcf_filter[5:-4]
                        try:
                            score_2 = vcf_record.genotype(name_2 + "---FIX---" + name_2_follow)["GT"]
                            group_score_2 +=float(score_2)
                            group_2_num +=1
                        except KeyError as ve_2:
                            error_2 = ve_2
                            pass
                try:
                    final_score = abs(group_score_1/group_1_num - group_score_2/group_2_num)
                except ZeroDivisionError:
                    final_score = 1
                    zero_record = "###No_cut_caused_group_missing_problem"
                    pass
                    

                if final_score > 0.75:
                    analysis_score += 1
        
        if not (error_1 == None):
            print("warning: " + str(error_1) + "cut with problem")
        if not (error_2 == None):
            print("warning: " + str(error_2) + "cut with problem")

        #print(vcf_filter[5:-4] + "---FIX---" + str(analysis_score))
        print(vcf_filter[5:-4] + "---FIX---" + str(analysis_score) + "---FIX---" + zero_record, file=open(args.OUT+".XIAO", "a"))
    for vcf_unfilter in vcf_unfilter_list:
        print(vcf_unfilter[5:-4] + "---FIX---NO_SNP---FIX---passed", file=open(args.OUT+".XIAO", "a") )  
    for vcf_unfilter_error in vcf_unfilter_error_list:
        print(vcf_unfilter_error[5:-4] + "---FIX---NO_CDS_GET---FIX---Unpassed", file=open(args.OUT+".XIAO", "a") )  
    subprocess.run(["rm", "-r", "./"+args.OUT], check=True)
    print("####################################################################")
    print("step 4.read and analysis the vcf files and get the scores passed")
    print("please check the file called " + args.OUT + ".XIAO")


if __name__ == '__main__':
    sys.exit(main())