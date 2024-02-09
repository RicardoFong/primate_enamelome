import re,sys,os,subprocess
from collections import defaultdict
from Functions_unmask import FASTA_iterator, translate

__author__= "Ricardo Fong Zazueta"
__mantainer__= "Ricardo Fong Zazueta"
__email__="ricardo.fong.zazueta@umontreal.ca"


'''
  This program is used to fix frameshifts in the problematic 
  sequences from the project.

  Usage:

  This_program.py Masked_file Trimmed_unmasked_fasta_file

  Where Masked_file is a tsv file with the following structure:

  Idividual_ID	Reference_genome	start-end,start-end

  Where start and end are the coordinates of nonsense variation that was manually masked.

  The program loops over the Masked_file and looks at the downstream and upstream of nonsense variation in the aligned but untrimmed files.
  when these coordinates are found, the exact coordinates of start and end are registered. It then looks for the corresponding nonsense
  sequence in the unaligned files and gets the start and end coordinates. Using this coordinates, the DNA file that yielded this protein
  sequence is edited by deleting one or two sequences in order to recover the frame of the desired protein.  
'''


# Directory with the unaligned fasta files 
Complete_unaligned=$UNALIGNED_SEQs_DIR

# Directory with the aligned files and trimmed
Trimmed=$TRIMMED_SEQs_DIR

# File with the coordinates that were masked manually
Masked_file=sys.argv[1]

# Trimmed_file=Trimmed+"/"+gene+"_trimmed_unmasked.mafft"
Trimmed_file=sys.argv[2]

# getting gene name
gene=Masked_file.split("_")[2]

with open(Masked_file, "r") as Masked_file:
	for line in Masked_file:
	
		data=FASTA_iterator(Trimmed_file)

		line=line.strip()
		ind=line.split("\t")[0]
		ref=line.split("\t")[1]

		# getting list of masked regions in sequence in a dictionary variable named coords
		coords=defaultdict(list)
		if "," in line.split("\t")[2]:

			thecoords=line.split("\t")[2].split(",")
			for i in thecoords:
				coord=i.split("-")
				coords[ind].append(coord)
		else:
			coords[ind].append(line.split("\t")[2].split("-"))

		tag=ind + "_" + gene + "." + ref # This tag is going to be used to find the sequence to be fixed within the trimmed unmasked files

		for seq in data:
			if tag in seq[0]: # If the sequence is in the masked files

				
				counter=0
				for coord in coords[ind]:
										
					Complete_fastas=Complete_unaligned + "/" + gene + "_complete.fasta"
					unaligned_data=FASTA_iterator(Complete_fastas)

					for ua_seq in unaligned_data:
						if tag in ua_seq[0]:

							header1=ua_seq[0]
							sequence1=ua_seq[1]

					# Here I am getting the last coordinate of the nonsense variation sequences. Exact, with no nucleotide variation
					strdata=FASTA_iterator(Trimmed_file)
					for strseq in strdata:
						if tag in strseq[0]:

							# Spotting the coordinates of the problematic sequences based on coords ###

							# This is the problematic nonsense sequence from the trimmed file, when there are gaps

							interest=seq[1][int(coord[0])-1:int(coord[1])]
							
							if "-" not in interest:

								# interest is the last 10 basepairs of the nonsense. you want to query the [1] coord
								interestend=interest[-10:]
								intereststart=interest[:10]
								
							# Searching the coordinates of the nonsense variation in the complete original sequence
							match=(re.search(interest, sequence1))

							if match is None:

								'''
								  Searching for coordinates of nonsense variation in original sequence
								  10 amino acids up and downstream are used, since the region can be trimmed, 
								  in which case, match can't be reached
								'''
								thematch=(re.search(interest,seq[1]))

								endinterest=thematch.span()[1]

																
								endinterest=seq[1][endinterest:endinterest+10]   # TWEAK, in case it is trimmed in this area
								endmatch=(re.search(endinterest, sequence1))
			
								if endmatch is None:
									
									endinterest=thematch.span()[1]
									endinterest=seq[1][endinterest-5:endinterest]
									endmatch=(re.search(endinterest, sequence1))
									
									# This match is used in case the interest region was trimmed and cannot be spotted
									end=endmatch.span()[1]
								else:
									endinterest=thematch.span()[1]
									endinterest=seq[1][endinterest:endinterest+10]

									endmatch=(re.search(endinterest, sequence1))
									
									# This is the normal match
									end=endmatch.span()[0]
								
								startinterest=thematch.span()[0]
								startinterest=seq[1][startinterest-8:startinterest-2] #TWEAK, in case it is trimmed in this area
								startinterest=startinterest.replace("-","")
								startmatch=(re.search(startinterest, sequence1))
								start=startmatch.span()[1]+2


								lastpartcoord=(end)*3
							else:
								
								interest=match.span()
								
								
								end=interest[1]			
								start=interest[0]

								# Here the end coordinate of the nonsense sequence is being spotted

								lastpartcoord=(end)*3

							# Building the first part of the sequence, all the way until the last aa before nonsense sequence
							
							if counter == 0:
								
								# Building the root sequence

								rootinterest=seq[1][:int(coord[counter])-1][-10:].replace("-","")

								match=(re.search(rootinterest, sequence1))
								
								therootinterest=match.span()
								
								# this is the root sequence before the first nonsense
								rootDNAcoord=therootinterest[1]*3

							else:
								'''
								   Building the midsequence
								   This is only for multiple problematic sites
								'''			
								midmatch=(re.search(prot_nonsense_seq,sequence1))
								midstart=midmatch.span()[1]
								
								probseq=seq[1][int(coord[0])-7:int(coord[0])-3]
								if "-" in probseq:													
									probseq=seq[1][int(coord[0])-7:int(coord[0])-1].replace("-","") # end of midsequence, located between two problematic sequecnces

								
								midendmatch=(re.search(probseq, sequence1))
								midend=midendmatch.span()[1]
								
								# Coordinates of the midsequence, protein and DNA 
								midcoord1=midstart
								midcoord2=midend

								DNAmidcoord1=midcoord1*3
								DNAmidcoord2=midcoord2*3
															
							break

					# This target is looking for the actual nonsense variation in the original sequence
					target=sequence1[start:end]
					
					target=re.search(target, sequence1).span()
					

					# This coordinate is exactly the beginning of the nonsense variation
					DNA_target=target[0]*3

					# Original DNA file		
					DNAfile=Complete_unaligned+"/"+gene+"_PGDP_DNA.fasta"
					

					DNAdata=FASTA_iterator(DNAfile)
					
					for DNAseq in DNAdata:
						if tag in DNAseq[0]:
							
							# First construct the root sequence:
							
							if counter==0:

								DNA_partseq=DNAseq[1][DNA_target:lastpartcoord]		   
								prot_nonsense_seq=translate(DNA_partseq)

								'''
								  add values to DNA_target and lastpartcoord in this ###
								  variable to edit the DNA sequence that yields the  ###
								  desired protein.
								'''
								DNA_partseq_final=DNAseq[1][DNA_target+1:lastpartcoord+1] # Tweak here


								prot_nonsense_seq_final=translate(DNA_partseq_final)
								
								rootDNA=DNAseq[1][:rootDNAcoord]
								rootseq=translate(rootDNA)

								DNA=rootDNA+DNA_partseq_final
								splitseq=rootseq + " " + prot_nonsense_seq_final
								
								if len(coords[ind]) == 1:

									DNArest=DNAseq[1][lastpartcoord:len(DNAseq[1])]
									
									DNA=rootDNA+ DNA_partseq_final + DNArest
									rest=translate(DNArest)

									splitseq=rootseq + " " + prot_nonsense_seq_final + " " + rest
									Protein=header1 + "\n" + translate(DNA)

							elif counter == 1:

								# Second construct the midsequence(s)
								DNA_partseq=DNAseq[1][DNA_target+1:lastpartcoord+1]			# Tweak here
								prot_nonsense_seq=translate(DNA_partseq)


								midDNA=DNAseq[1][DNAmidcoord1:DNAmidcoord2]
								midProt=translate(midDNA)
							
								DNArest=DNAseq[1][lastpartcoord:len(DNAseq[1])]
							
								rest=translate(DNArest)

								splitseq = header1+"\n"+splitseq + " " + midProt + " " + prot_nonsense_seq + " " + rest
									
								DNA=DNA + midDNA + DNA_partseq + DNArest
								
								Protein=header1+"\n"+translate(DNA)
					
							counter+=1

							break

				newDNAfile=$DNA_OUT_FILE
				newProtfile=$PROT_OUT_FILE
				

				# Write edited DNA and Protein sequences to files

				with open(newDNAfile, "a") as DNAfile:
					DNAfile.write(DNAseq[0] + "\n" + DNA+"\n")
				
				with open(newProtfile, "a") as Protfile:
					Protfile.write(Protein +"\n")

				print(splitseq+"\n")
				print(Protein+"\n")
				break
