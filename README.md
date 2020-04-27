# Corevar
bacteria core-genome comparison betweenn different clades.
# introduction
Corevar is a  is a workflow for core-genome comparison betweenn different clades from raw sequencing reads or genome assemblies.
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
# Quick Start (run it on Linuc system)
## xiao_robot_match_classify_CDS.py

## xiao_robot_extract_CDS.py

## xiao_extract_beside_CDS.py

## xiao_robot_SNP_analysis_between_groups.py

