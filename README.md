# NARASIMHA Primer Designer

## Introduction
A program to generate primers for the NARASIMHA Assay for Targeted RNA Sequencing to identify Chimeric Gene Fusions: https://www.nature.com/articles/s41408-020-0313-6  

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
   Create a 6 column BED file with the follwing values in tab delimited columns.
       chromosome,  start-coordinate  end-coordinate,  unique gene/exon name (without spaces),  Score(default 1 for all),   strand(+ or -)
   Create separate BED files for plus and minus strand genes
       
       EXAMPLE
         chr9  133710452  133710912  ABL1_EX1  1  +
         chr9  133729450  133729624  ABL1_EX2  1  +

~OR~
~###    b. FASTA/multi-FASTA file~
   ~Alternatively, a fasta file can be provided as input with cDNA sequence of the gene of interest.~


##  3. RUN
~python3 ./script/narasimha.py -i <path to BED or FASTA input file> -g <path to reference genome (.fasta)>~

    python3 ./script/narasimha_plus.py -i <path to BED file with + strand genes> -g <path to reference genome (.fasta)> -l <length of region for primer design> 
    python3 ./script/narasimha_minus.py -i <path to BED file with - strand genes> -g <path to reference genome (.fasta)> -l <length of region for primer design>
    EXAMPLE:
    ./narasimha_plus.py -i lung_panel_v3_format_plus.bed -g /home/reference_genomes/hg19_ref_index/hg19.fasta -l 65
    ./narasimha_minus.py -i lung_panel_v3_format_minus.bed -g /home/reference_genomes/hg19_ref_index/hg19.fasta -l 65

## 4. Output
    A "primer_out_*" directory will be created as output. The narasimha_minus.py creates primer_out_minus as output directory, narasimha_plus.py creates primer_out_plus as output directory.
    The round1 and round2 primers are present in the "final" directory present inside the above folders.
    The final list of primers for both rounds can be found in "all_primers.csv" in the "final" directory.
    If no primers are generated, this may be due to short length of sequences.  
    ~For sequences less than 50 bases, it is recomended that the sequence of the previous exon also be taken and a fasta file be created. Run the program with the fasta file as input.~

## 5. Visualization
  The scripts 


