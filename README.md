# NARASIMHA Primer Designer

## Introduction
A program to generate primers for the NARASIMHA Assay: https://www.nature.com/articles/s41408-020-0313-6

## 1. Prerequisites
python=3.6  
pandas=1.0.1  
biopython  
bedtools=2.29.1  
primer3 (libprimer3 release 2.5.0)  
isPcr

or create and activate environment in scripts directory  
```
        $ conda env update -f ./scripts/environment.yaml  
        $ conda activate narasimha-primers
```


##  2.  Creating Input file
###    a. BED file
   Create a 6 column BED file with the follwing values in tab delimited columns
       
       chromosome,  start-coordinate  end-coordinate,  unique gene/exon name (without spaces),  Score(default 1 for all),   strand(+ or -)
       
       EXAMPLE
         chr9  133710452  133710912  ABL1_EX1  1  +
         chr9  133729450  133729624  ABL1_EX2  1  +

OR
###    b. FASTA/multi-FASTA file  
   Alternatively, a fasta file can be provided as input with cDNA sequence of the gene of interest.


###  3. RUN

    python3 ./script/narasimha.py -i <path to BED or FASTA input file> -g <path to reference genome (.fasta)>


### 4. Output
    A "FINAL" directory will be created. The final list of primers for both rounds can be found in "all_primers.csv" in this directory.  
    If no primers are generated, this may be due to short length of sequences.  
    For sequences less than 50 bases, it is recomended that the sequence of the previous exon also be taken and a fasta file be created. Run the program with the fasta file as input.  





