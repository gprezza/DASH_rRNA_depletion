# DASH_rRNA_depletion
Python tool to design sgRNA pools targeting the rRNA genes, as described in [this paper](https://rnajournal.cshlp.org/content/26/8/1069.long) (PMID: 32345633). 

Example input and output files for *Salmonella* Typhimurium (gff annotation) and *Bacteroides thetaiotaomicron* (custom rRNA annotation with the "bt_rRNA.bed" file) are present in the "example_files" folder.

## Dependencies
Create a conda environment from the environment.yml file:
```
conda env create -f environment.yml
```
This will install all required dependencies in an environment called "DASH".

## Usage
```
usage: design.py [-h] [--gff GFF] [--manual_ann MANUAL_ANN] [--minGC MINGC]
                 [--maxGC MAXGC] [--length LENGTH] [--offtargets] [--pam PAM]
                 fasta_file

Designs a pool of sgRNAs targeting the ribosomal genes of a bacterial species.
Indicate an annotation file with either the --gff or --manual_ann options

positional arguments:
  fasta_file            Fasta file with the genome sequence

optional arguments:
  -h, --help            show this help message and exit
  --gff GFF, -g GFF     GFF file with the genome annotation (Default: False).
  --manual_ann MANUAL_ANN, -ma MANUAL_ANN
                        Use if you want to provide a file with the manual
                        annotation (in bed format, tab-separated) of the rRNA
                        genes (Default: False).
  --minGC MINGC, -gc MINGC
                        Minimal accepted GC% of a spacer (Default: 30).
  --maxGC MAXGC, -GC MAXGC
                        Maximal accepted GC% of a spacer (Default: 80).
  --length LENGTH, -l LENGTH
                        Spacer length (Default: 20).
  --offtargets, -o      Print the spacers that were discarded because of off-
                        targeting.
  --pam PAM, -p PAM     PAM sequence (Default: NGG).
```
### Input

The working directory must have the **design_grnas.py** file.

* If you already have Bowtie indexes of these sequences, put them in a subfolder named bowtie_files/. If you don't, the subfolder and the indexes will be generated by the script.

* If you want to use a gff3 annotation, give the path after the --gff option. The rRNA genes must have the "rRNA" type (3rd column of the gff file). The scaffold ID(s) in column 1 of the file must be the same present in the genome fasta file.

* If you don't want to use a gff annotation, provide the coordinates of the rRNA genes **in BED format** through the --manual_ann argument. The file must have the following fields (tab-separated and in this order; do not include a header): scaffold, start, end, name, score, strand
  - scaffold: ID of the scaffold (chromosome or plasmid). Must be the same present in the respective fasta file.
  - start: start position of the rRNA gene.
  - end: end position of the rRNA gene.
  - name: name of the rRNA gene. Must contain one among: 5S, 16S or 23S.
  - score: required by the BED format but not used by the script. Can be anything (e.g. a dot).
  - strand: strand of the rRNA gene.

  An example file is the following, for *B. thetaiotaomicron* rRNAs:
  ```
  NC_004663.1	1626501	1626612	5S_r01	.	-
  NC_004663.1	1626629	1629600	23S_r02	.	-
  NC_004663.1	1630071	1631660	16S_r03	.	-
  NC_004663.1	2336964	2337075	5S_r04	.	-
  NC_004663.1	2337153	2340121	23S_r05	.	-
  NC_004663.1	2340595	2342191	16S_r06	.	-
  NC_004663.1	3030969	3031080	5S_r07	.	-
  NC_004663.1	3031158	3034127	23S_r08	.	-
  NC_004663.1	3034601	3036195	16S_r09	.	-
  NC_004663.1	3325219	3325330	5S_r10	.	-
  NC_004663.1	3325408	3328377	23S_r11	.	-
  NC_004663.1	3328851	3330446	16S_r12	.	-
  NC_004663.1	4983424	4985021	16S_r13	.	+
  NC_004663.1	4985493	4988462	23S_r14	.	+
  NC_004663.1	4988479	4988590	5S_r15	.	+
  ```

After all files have been copied, run the **design_grnas.py** script.

**NB** To maximize depletion efficiency, it's important that the rRNA genes are annotated as precisely as possible. If you already have RNA-seq data of your species of interest, it's a good idea to check if the extremities of the rRNA genes are annotated correctly in the gff3/bed file and modify the file accordingly. 

### Output
The following files are generated by the script:

* oligos.csv: file with the sequences of the oligos that have to be ordered to generate the pool.
* grnas.fa: multifasta file with all the spacers selected by the script.
* bowtie.csv: output file of Bowtie (to have a look at oligos discarded for off-targeting).

