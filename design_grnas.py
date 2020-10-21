from Bio import Seq
from Bio import SeqIO
from itertools import product
import sys
import csv
import subprocess
import re
import os
import argparse
import primer3 #primer3-py package

parser = argparse.ArgumentParser(description='Designs a pool of sgRNAs '
								'targeting the ribosomal genes of a bacterial species.')

parser.add_argument('--minGC', '-g', default=30, type=int,
					help='Minimal accepted GC%% of a spacer (Default: 30).')

parser.add_argument('--maxGC', '-G', default=80, type=int,
					help='Maximal accepted GC%% of a spacer (Default: 80).')

parser.add_argument('--length', '-l', default=20, type=int,
					help='Spacer length (Default: 20).')

parser.add_argument('--offtargets', '-o', default=False, action='store_true',
					help='Print the spacers that were discarded because of off-targeting.')

parser.add_argument('--manual_ann', '-ma', default=False, action='store',
					help='Use if you want to provide a file with the manual '
					'annotation (in bed format) of the rRNA genes (Default: False).')

parser.add_argument('--pam', '-p', default="NGG", type=str,
					help='PAM sequence (Default: NGG).')


args = parser.parse_args()
GC_low = args.minGC
GC_high = args.maxGC
guide_length = args.length
show_offtargets = args.offtargets
annotation_file = args.manual_ann
PAM = args.pam.upper()

ambiguous_alph = {"N":["A","C","G","T"], "V":["A","C","G"], "H":["A","C","T"], "D":["A","G","T"],
	"B":["C","G","T"], "M":["A","C"], "K":["G","T"], "W":["A","T"], "S":["C","G"], "Y":["C","T"],
	"R":["A","G"], "A": ["A"], "C": ["C"], "G": ["G"], "T": ["T"]}

def remove_ambiguous(sequence):
	#returns a list of all possible sequences
	#given the ambiguous sequence input (e.g. with N,V, etc. bases)
	return ["".join(i) for i in product(*[ ambiguous_alph[j] for j in sequence ]) ]

def revcomp(sequence):
#returns reverse complement of sequence
	complement = {'A':'T', 'C':'G', 'T':'A', 'G':'C', 'N':'N','Y':'R',
	'R':'Y', 'W':'W', 'S':'S', 'K':'M', 'M':'K', 'D':'H', 'H':'D',
	'V':'B', 'B':'V', 'X':'X', '-':'-'}
	string_list = list(sequence)
	string_list.reverse()
	revcompl = ''
	for base in string_list:
		revcompl += complement[base]
	return revcompl

def find_all(string, substr):
#makes a list with the indexes of all occurrencies of substr in string
	listindex=[]
	start = 0
	i = string.find(substr, start)
	while i >= 0:
		listindex.append(i)
		i = string.find(substr, i + 1)
	return listindex

def GC_count(string):
#returns %GC of string
	bases = list(string)
	G = bases.count("G")
	C = bases.count("C")
	GC_content = 100*(G+C)/len(string)
	return GC_content

def retrieve_guide(PAM,strand):
#get spacer sequence given the index of the pam sequence (PAM)
#and the strand of the pam
	if strand == "bot":
		start = PAM - guide_length
		if start > 0:
			end = PAM
			guide = sequence[start:end]
		else: #if the spacer would start outside the rRNA gene
			guide = None
	else:
		start = PAM + PAM_length
		end = start + guide_length
		if end < len(sequence):
			guide = revcomp(sequence[start:end])
		else: #if the spacer would end outside the rRNA gene
			guide = None
	return guide

def offtarget(list_of_PAMs):
#checks if the offtarget is followed by a PAM
#and if it's outside a rRNA gene.
#If so, removes the spacer from guides_dict and returns 1.
	if any (genome_seqs[scaffold][start:end] == PAM for PAM in list_of_PAMs):
		if scaffold in positions:
			if not any(index in range(i,j+1) for [i, j] in positions[scaffold]):
			#if the index is not in at least one rRNA gene:
				guides_dict.pop(row[0], None)
				n = 1
				if show_offtargets == True:
					print("\t".join(row))
			else:
				n = 0
		else:
			n = 0
	else:
		n = 0
	return n

def primer_heterodimer(oligo,ID):
#checks for unwanted internal binding of the fill-in primer to oligo.
#cutoffs are predicted melting temperature of the annealing >40 °C and 
#of the 3' portion of 30 °C
	res_end = primer3.bindings.calcEndStability(primer_revcomp,oligo)
	res_het = primer3.bindings.calcHeterodimer(primer_revcomp,oligo,output_structure=True)
	if res_het.tm > 40 or res_end.tm > 30:
		print("\nThe end-filling primer has high predicted internal binding to the", ID, "oligo\nThis oligo is discarded.")
		print("Heterodimer formation:")
		print(res_het)
		print("End Stability:")
		print(res_end)
		print("Predicted binding:")
		print(res_het.ascii_structure)
		ID = None
	return ID

#cleanup of previous interrupted runs:
if os.path.isfile('bowtie.csv'):
        os.remove('bowtie.csv')
if os.path.isfile('grnas.fa'):
        os.remove('grnas.fa')
if os.path.isfile('reference_sequences/genome_sequence.fa'):
        os.remove('reference_sequences/genome_sequence.fa')
if os.path.isfile('input.fa'):
        os.remove('input.fa')

#check if all required files exist and create folders if needed:
if annotation_file is False:
	if len(os.listdir("annotations")) == 0:
		print("No annotation file in the annotations folder! Did you mean "
		"to use the --manual_ann argument?")
		sys.exit()

if not os.path.isdir('bowtie_files'):
	os.mkdir('bowtie_files')


if len(os.listdir("reference_sequences")) == 0:
	print("No files are present in the reference_sequences "
	"folder, exiting the script.")
	sys.exit()

#make a multifasta file with all sequences present in /reference_sequences
genome_fa_file = open("reference_sequences/genome_sequence.fa",mode= "w+")
for filename in os.listdir("reference_sequences"):
	if filename != "genome_sequence.fa":
		_file = open("reference_sequences/{0}".format(filename),mode= "r+")
		for line in _file:
			genome_fa_file.write(line)
		genome_fa_file.write("\n")
genome_fa_file.seek(0)

#make a unique gff file with all files present in /annotations
if annotation_file is False:
	gff_file = open("annotations/annotation_file.gff",mode= "w+")
	for filename in os.listdir("annotations"):
		if filename != "annotation_file.gff":
			_file = open("annotations/{0}".format(filename),mode= "r+")
			for line in _file:
				gff_file.write(line)
	gff_file.seek(0)

regex_rrna = re.compile('\d*S') #matches rRNA name (5S, 16S, 23S)

PAM_length = len(PAM)

T7 = "TTCTAATACGACTCACTATA" #T7 promoter (minus the first G)
scaff = "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTT"
#scaff is the Cas9 sgRNA scaffold
primer = "GGCATACTCTGCGACATCGT" #primer for fill-in of the oligos. The reverse complement
				#gets added to the sgRNA template oligos.

###########main:

total = 0 #no. of grnas that passed all criteria
one_G = 0 #spacers starting with one G
two_G = 0 #spacers starting with two G
total_nt = 0 #total no. of nucleotides in the final oligo pool

PAM_list = remove_ambiguous(PAM)

#create list of reverse complement of PAMs (for off-target analysis later)
revcomp_PAM_list = []
for PAM in PAM_list:
	PAM = revcomp(PAM)
	revcomp_PAM_list.append(PAM)

rRNA_genes = {"5S":[], "16S":[], "23S":[]}
positions = {} #dictionary of start end end positions of the rRNA genes
				#scaffold:[[start1,end1],[start2,end2]]

print("Looking for spacers.")

#create genome_seqs dictionary with sequences of chromosome and plasmids,
#"positions" dictionary with coordinates of all rRNA genes
#and rRNA_genes dictionary with sequences of all rRNA genes:
genome = SeqIO.parse(genome_fa_file, 'fasta')
genome_seqs = {}
for record in genome:
	if "genome" in record.description: #if it's the main chromosome
		genome_id = re.findall('[^ ]+', record.description)[0]
		genome_seqs[genome_id] = str(record.seq)
	else:
		plasmid_id = re.findall('[^ ]+', record.description)[0]
		genome_seqs["{0}".format(plasmid_id)] = str(record.seq)
		#adds the plasmids to the genome_seqs dictionary

if annotation_file is False: #if a gff file is in the \annotations folder
	for line in gff_file:
		if line[0] != "#":
			linelist = line.split("\t")
			if linelist[2] == "rRNA":
				ID = regex_rrna.findall(linelist[8])[0]
				start = int(linelist[3])
				end = int(linelist[4])
				strand = linelist[6]
				scaffold = linelist[0]
				positions.setdefault(scaffold,[]).append([start,end])
				if strand == "-":
					seq = revcomp(genome_seqs[scaffold][start:end])
				else:
					seq = str(genome_seqs[scaffold][start:end])
				rRNA_genes[ID].append(seq)
else: #if a custom rRNA annotation was provided
		annotation = open("{0}".format(annotation_file), "r")
		reader = csv.reader(annotation, delimiter = "\t")
		for line in reader:
			try:
				scaffold = line[0]
				ID = regex_rrna.findall(line[3])[0]
				start = int(line[1])
				end = int(line[2])
				positions.setdefault(scaffold,[]).append([start,end])
				#line above appends to key if already exists, otherwise creates
				#key with [[start, end]] as value
				strand = line[5]
				if strand == "-":
					seq = revcomp(genome_seqs[scaffold][start:end])
					rRNA_genes[ID].append(seq)
				else:
					seq = genome_seqs[scaffold][start:end]
				rRNA_genes[ID].append(seq)
			except:
				print("Error: something went wrong while handling the custom "
				"rRNA annotation file provided after the -ma argument. Is the "
				"file in the right (bam) format?")
				sys.exit()

#find gRNA spacer sequences in the rRNA genes and filter for GC content:
invalid_GC = [] #list of guides discarded because of GC content
guides = [] #list of lists. Each sublist is [ guide sequence, rRNA gene targeted ]
for rRNA in rRNA_genes:
	for sequence in rRNA_genes[rRNA]: #loops through the sequence of each copy of all rRNAs
		PAMs_top = [] #PAM sequences in the top strand
		PAMs_bot =  [] #PAM sequences in the bottom strand
		for PAM_sequence in PAM_list:
			PAMs_bot.extend(find_all(sequence, PAM_sequence))
			PAMs_top.extend(find_all(sequence, revcomp(PAM_sequence)))

		for PAM in PAMs_bot:
			guide = retrieve_guide(PAM, "bot")
			if guide and [guide, rRNA] not in guides:
				if GC_low <= GC_count(guide) <= GC_high:
					guides.append([guide, rRNA])
				else:
					if [guide, rRNA] not in invalid_GC:
						invalid_GC.append([guide, rRNA])
		for PAM in PAMs_top:
			guide = retrieve_guide(PAM, "top")
			if guide and [guide, rRNA] not in guides:
				if GC_low <= GC_count(guide) <= GC_high:
					guides.append([guide, rRNA])
				else:
					if [guide, rRNA] not in invalid_GC:
						invalid_GC.append([guide, rRNA])


#search for off-targets and discard those gRNAs:
print("Looking for off-targets.")
guides_dict = {} #dictionary of gRNA_ID : sequence
i = 1
previous_rrna = ""
with open('input.fa', 'w+') as f: #input file for bowtie
	for guide in guides:
		if guide[1] == previous_rrna:
			i += 1
		else:
			i = 1
		f.write(">")
		line = guide[1] + "_" + str(i) #i is the count of gRNAs targeting
						#the same rRNA gene
		f.write(line) #ID of the fasta entry: targeted rRNA + _i
		f.write("\n")
		f.write(guide[0]) #gRNA sequence
		f.write("\n")
		guides_dict[line] = guide[0]
		previous_rrna = guide[1]


#make bowtie indexes if not present:
if not os.path.isfile('bowtie_files/*.ebwt'):
	cmd = ("bowtie-build -r -q -f reference_sequences/genome_sequence.fa bowtie_files/genome_index"),
	subprocess.call(cmd, shell=True, stderr=subprocess.DEVNULL)

#run bowtie:
offtargeting = 0 #no. of gRNAs discarded because of off-targeting
with open("bowtie.csv", mode='w+') as bowtie_out: #output of bowtie with all off-targets
	cmd = ("bowtie -v 3 -a bowtie_files/genome_index --suppress 6,7 "
	"-f input.fa >> bowtie.csv") #matches with up to 3 mismatches
	subprocess.call(cmd, shell=True, stderr=subprocess.DEVNULL)
	csv_reader = csv.reader(bowtie_out, delimiter='\t')
	if show_offtargets == True:
		print("\nFound and discarded these offtargets:\n")
		print("name\tstrand\tscaffold\tstart index\tsequence\tmismatches")
	for row in csv_reader:
		scaffold = row[2] #chromosome or plasmid ID
		index = int(row[3]) #starting index of off-target
		if row[1] == "+": #strand
			start = index + guide_length
			end = index + guide_length + PAM_length
			offtargeting += offtarget(PAM_list) #if offtarget() identifies
							#this sequence as a true
							#off-target, the gRNA is
							#discarded and offtargeting
							#is increased by 1. Otherwise
							#by 0.
		else:
			start = index - PAM_length
			end = index
			offtargeting += offtarget(revcomp_PAM_list)

#make the oligonucleotide sequences and check for unwanted binding of fill-in primer:
one_G = 0 #spacers starting with one G
two_G = 0 #spacers starting with two G
total_nt = 0 #total no. of nucleotides in the final oligo pool
primer_revcomp = revcomp(primer)

print("Looking for unwanted internal binding of the fill-in primer to each oligo.")
with open("grnas.fa","w+") as grnas_out: #file with all final gRNAs
	with open("oligos.csv", "w+") as output: #file with all oligonucleotides
		writer = csv.writer(output)
		for gRNA_id, sequence in guides_dict.items():
			if sequence[0] == "G":
				if sequence[1] == "G": #starts with two G
					oligo = T7 + sequence + scaff + primer_revcomp
					ID = primer_heterodimer(oligo[:-len(primer_revcomp)],gRNA_id)
					if ID is not None:
						grnas_out.write(">")
						grnas_out.write(gRNA_id)
						grnas_out.write("\n")
						grnas_out.write(sequence)
						grnas_out.write("\n")
						line = [gRNA_id, oligo]
						two_G += 1
				else: #starts with only one G
					oligo = T7 + "G" + sequence + scaff + primer_revcomp
					ID = primer_heterodimer(oligo[:-len(primer_revcomp)],gRNA_id)
					if ID is not None:
						grnas_out.write(">")
						grnas_out.write(gRNA_id)
						grnas_out.write("\n")
						grnas_out.write(sequence)
						grnas_out.write("\n")
						line = [gRNA_id, oligo]
						one_G += 1
			else:
				oligo = T7 + "GG" + sequence + scaff + primer_revcomp
				ID = primer_heterodimer(oligo[:-len(primer_revcomp)],gRNA_id)
				if ID is not None:
					grnas_out.write(">")
					grnas_out.write(gRNA_id)
					grnas_out.write("\n")
					grnas_out.write(sequence)
					grnas_out.write("\n")
					line = [gRNA_id, oligo]
			writer.writerow(line)
			total_nt += len(oligo)

total = len(guides_dict) #no. of grnas that passed all criteria
zero_G = total - one_G - two_G
print("Done! Some (maybe) useful data:\n\n"
"The final pool consists of {0} gRNAs.\n\n"
"{1} gRNAs were discarded on the basis of GC content.\n"
"{2} gRNAs were discarded for off-targeting.\n"
"{3} gRNAs don't start with G.\n"
"{4} gRNAs start with one G.\n"
"{5} gRNAs start with two G.\n"
"There is a total of {6} nucleotides."
"".format(total, len(invalid_GC), offtargeting, zero_G, one_G, two_G, total_nt))

#cleanup
os.remove("reference_sequences/genome_sequence.fa")
if annotation_file is False:
	os.remove("annotations/annotation_file.gff")
os.remove("input.fa")
