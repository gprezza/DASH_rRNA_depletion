# DASH_rRNA_depletion
Python tool to design sgRNA pools targeting the rRNA genes.

## Dependencies
These Python (3.x) packages are required to run th the script:

* Bio
* primer3-py

You'll also need to have [Bowtie](http://bowtie-bio.sourceforge.net/index.shtml) installed.

## Usage
```
Designs a pool of sgRNAs targeting the ribosomal genes of a bacterial species.

optional arguments:
  -h, --help            show this help message and exit
  --minGC MINGC, -g MINGC
                        Minimal accepted GC% of a spacer (Default: 30).
  --maxGC MAXGC, -G MAXGC
                        Maximal accepted GC% of a spacer (Default: 80).
  --length LENGTH, -l LENGTH
                        Spacer length (Default: 20).
  --offtargets, -o      Print the spacers that were discarded because of off-
                        targeting.
  --manual_ann MANUAL_ANN, -ma MANUAL_ANN
                        Use if you want to provide a file with the manual
                        annotation (in bed format) of the rRNA genes (Default:
                        False).
  --pam PAM, -p PAM     PAM sequence (Default: NGG).
```
### Input

The working directory must have the **design_grnas.py** file and the subfolders **reference_sequences/** and **bowtie_files/**. If you want to use a gff3 annotation of the transcriptome, you'll also need the **annotations/** folder.

* Copy the fasta sequences of the chromosome and plasmids (if any) in the reference_sequences/ folder. 

* If you already have bowtie indexes of these sequences, copy them into bowtie_files/, otherwise they will be generated by the script.

* If you want to use a gff3 annotation, copy the file(s) in the annotations/ folder. The rRNA genes must have the "rRNA" type (3rd column of the gff file). The scaffold ID(s) in column 1 of the file must be the same present in the fasta file(s) of reference_sequences/.

* If you don't want to use a gff annotation, provide the coordinates of the rRNA genes **in BED format** through the -manual_ann argument. The file must have the following fields (tab-separated and in this order; do not include a header): scaffold, start, end, name, score, strand
  - scaffold: ID of the scaffold (chromosome or plasmid). Must be the same present in the respective fasta file.
  - start: start position of the rRNA gene.
  - end: end position of the rRNA gene.
  - name: name of the rRNA gene. Must contain one among: 5S, 16S or 23S.
  - score: required by the BED format but not used by the script. Can be anything (e.g. a dot).
  - strand: strand of the rRNA gene.

  An example file is the following, for _B. thetaiotaomicron_ rRNAs:
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

