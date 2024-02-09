from Functions import FASTA_iterator
import sys, os

__author__= "Ricardo Fong Zazueta"
__mantainer__= "Ricardo Fong Zazueta"
__email__="ricardo.fong.zazueta@umontreal.ca"

'''
	This program is used to fragment the complete protein sequences based on experimental data.
	The coordinates used to fragment this data were obtained from alignments of peptide sequences against 
	the computationally obtained sequences per protein.

	Functioning:
	python3 thisprogram.py

	The directories in the header section must be set for proper functioning with the subset coordinate files and fasta files
	per protein in-place

	Simulated Fragmented sequences per proteins are written to $outdir


'''

subsetDir=$subDir # directory with the coordinates of experimental data

fastaDir=$fDir # directory with fasta alignments per protein

outdir=$oDir # output directory

sub_files=os.listdir(subsetDir)

genes=[]
ages=[]


# Building gene and ages lists from the experimental coordinates files
for file in sub_files:
	thyseq=[]
	gene=file.split("_")[0]
	age=file.split("_")[1].split(".")[0]
	ageFile=subsetDir+"/"+file
	
	protein=fastaDir+"/"+gene+".fasta" # aligned fasta file of single proteins
	
	data=FASTA_iterator(protein)
	outfile=gene+"_"+age+".fasta"

	thyseq=[]
	
	print("Writting subset sequences of gene "+ gene + " of temporality " + age +" to file:\n" + outfile)
	
	for seq in data:
		
		data=FASTA_iterator(protein)
		thyseq=[]	
		with open(ageFile, "r") as aFile:
			for line in aFile:
				line=line.strip()
				
				start=int(line.split("-")[0])
				end=int(line.split("-")[1])
				partseq=seq[1][start:end] # subsetting sequence by start-end coordinate

				thyseq.append(partseq)

		outfile=outdir+"/"+gene+"_"+age+".fasta" # Naming the outfile
		
		subsetseq="".join(thyseq) # This sequence is constructed by all partial sequences ($partseq) of each individual coordinate

		# writting fragmented sequence by protein file
		with open(outfile, "a") as thyoutfile:

			thyoutfile.write(seq[0] + "\n" + subsetseq + "\n")
		
