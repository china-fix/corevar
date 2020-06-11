# Corevar
bacteria core-genome comparison between different clades.
# Introduction
Corevar is a workflow for core-genome comparison between different clades from raw sequencing reads or genome assemblies.
# Dependencies
1. [python 3.6.9](https://www.python.org/) and following librarires are required:
* [Biopython 1.74](https://biopython.org/)
* [PyVCF 0.6.5](https://pyvcf.readthedocs.io/en/latest/index.html)
2. [SPAdes 3.12.0](http://cab.spbu.ru/software/spades/);
3. [Prokka 1.13.3](https://github.com/tseemann/prokka);
4. [Roary 3.12.0](https://sanger-pathogens.github.io/Roary/);
5. [MAFFT 7.310](https://mafft.cbrc.jp/alignment/software/);
6. [PRANK 170427](http://wasabiapp.org/software/prank/);
7. [SNP-sites](https://github.com/sanger-pathogens/snp-sites);
8. [blastn 2.6.0+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
# Corevar workflow
![](workflow.jpg?raw=true "Corevar workflow")
* Flowcharts of the steps in the application of the Corevar. (A) Unique CDS analysis workflow. (B) Core-SNP and core-upstream-SNP analysis workflow.
# Quick Start (run it on Linux system)
## xiao_robot_match_classify_CDS.py
* This is a python script to locally blastn the CDS.fasta to the reference genome and filter the orginal CDS.fasta to matched_CDS.fasta and un_matched_CDS.fasta  

Usage: python3 xiao_robot_match_classify_CDS.py

    --CDS (required=True, type=str, metavar='FILENAME', help="the CDS fasta filename you want to blast")
    
    --REF (required=True, type=str, metavar='FILENAME', help="the complete reference genome of fasta format you want to blast to")
    
    --CUTOFF (default=0.9, type=float, metavar='DEFAULT 0.9', help="the lowest similarity value which classify as matched")
    
    --OUT (default="xiao_robot_match_classify_CDS", type=str, metavar='directory', help="Output directory name")

## xiao_robot_extract_CDS.py
* This python script is modified from the script of xiao_robot_extract_beside_CDS. Script do blastn and extract the CDS in each strain.

Usage: python3 xiao_robot_extract_CDS.py

    --CDS (required=True, type=str, metavar='FILENAME', help="the CDS fasta filename you want to extract")
    
    --DB (required=True, type=str, metavar='FILENAME', help="the database name in the command blastn -db")
    
    --WASH_LIST (required=True, type=str, metavar='FILENAME', help="the filename (bacteria strain fasta assemblies) list you want to get the CDS")
    
    --OUT (default="xiao_robot_extract_CDS", type=str, metavar='directory', help="Output directory name")

## xiao_extract_beside_CDS.py
* This python script is used to blastn the input CDS fasta file to the assembly database, then get the location and other hsp information. Finally, extract the upstream CDS sequences in every bascteria assemblies.

Usage: pyhon3 xiao_extract_beside_CDS.py

    --CDS (required=True, type=str, metavar='FILENAME', help="the CDS fasta filename you want to blast")
    
    --DB (required=True, type=str, metavar='FILENAME', help="the database name in the command blastn -db")
    
    --WASH_LIST (required=True, type=str, metavar='FILENAME', help="the filename (bacteria strain fasta assemblies) list you want to get the upstream sequence of relative CDS")
    
    --CUT (default=100, type=int, metavar='up-stream cut length', help="the length(bp) you want to cut upstream of the CDS")
    
    --OUT (default="xiao_robot_extract_beside_CDS", type=str, metavar='directory', help="Output directory name")

## xiao_robot_SNP_analysis_between_groups.py
* This is a python script used to analysis the SNP between different bacteria strain clades. Trying to find the conserved SNP (core-SNP) cases between the user defined clades. 
1. extract core CDS in every strains according core CDS reference
2. multiple sequence alignment of each core CDS for all the strains by the program PRANK or MAFFT
3. SNP variation analysed by SNP-sites
4. calculate score

Usage: python3 xiao_robot_SNP_analysis_between_groups.py

    --CDS_LIST (required=True, type=str, metavar='FILENAME', help="the CDS list you want to extract")
    
    --STRAIN_LIST (required=True, type=str, metavar='FILENAME', help="the strain list you want to extract from the fasta files (each file comtain all the CDS of a strain)")
    
    --GROUP_1 (required=True, type=str, metavar='FILENAME', help="the strain (fasta files name) list you want to caculate as one clade")
    
    --GROUP_2 (required=True, type=str, metavar='FILENAME', help="the strain (fasta files name) list you want to caculate as one clade")
    
    --TEMP_SAVE (action='store_const', const=True, metavar='SAVE_SOME_TEMP', help="this command help to save some temp infomation or files")
    
    --FAST (action='store_const', const=True, metavar='FAST ALIGNMENT WITH MAFFT', help="this command help to use MAFFT to do a fast alignment")
   
    --OUT (default="xiao_robot_SNP_analysis_between_groups", type=str, metavar='directory', help="Output directory name")
# Citation

 
