** This code is explicitly for conversion and extraction of protein-coding genes when aligned to human genome assembly build 38 v23 using STAR. **

** For those users using a different human genome build, please make sure that your input genes are in the ensembl format. We have provided a .csv file with ensemble ID and gene names called "ENSG_GeneSymbol_list.csv" ('ensembl' is used as the final input for the classifier) . **


This folder contains a python script, extract_PCG.py to take in a tab separated .txt file of gene counts from one or more samples (ReadsPerGene.out.tab, exported from STAR) and export only the protein coding genes that is used as the input for the classifier. 


Demo
Usage: python extract_PCG.py input_file.txt output_file.csv

The output file is a .csv file with 19,590 protein-coding genes used as input for MLC1 and MLC2. The expected run time for demo on normal desktop computer is 5.863s.


Instructions for use

python extract_PCG.py data_files/test1_combined_data.txt data_files/test_set_1.csv

The data_files folder contains the mock_combined_data.txt which is a tab delimited file with gene count information across 114 samples. This python script uses ENSG_GeneSymbol_list.csv, also provided in the same folder, to identify the protein coding genes from the inputted .txt file.

