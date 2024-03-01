from Functions import FASTA_iterator
import sys, os

__author__= "Ricardo Fong Zazueta"
__mantainer__= "Ricardo Fong Zazueta"
__email__="ricardo.fong.zazueta@umontreal.ca"

'''

	This program is used to build the subset files of rate for site, and shannon entropy over, and under the threshold of 1.
	The folers where these subset files are built and stored are the ones named:

	r4sDir
	shannonDir

	The functioning of the program is:
	python3 thisprogram.py $num

	where $num is a string value, any of 5, 10 or 14, refering to the number of proteins in each file concatenation.

	The program outputs:
	1.- Fasta files of the concatenation defined in $num, subset by Rate4site values over, and under threshold of 1
	2.- Fasta files of the concatenation defined in $num, subset by Shannon entropy values over, and under threshold of 1
	3.- Nexus files with the coordinates of aminoacids belonging to each protein to be used in IQtree for over and under threshold of 1 for Rate4site
	4.- Nexus files with the coordinates of aminoacids belonging to each protein to be used in IQtree for over and under threshold of 1 for Shannon entropy

'''

# Declaring input files

ProtNum=sys.argv[1] # one number related to concatenations (5, 10 or 14)

r4sDir=$r4sDir # directory with the Rate4site values 
shannonDir=$shDir # directory with the shannon entropy values
basedir=$fastaDir # directory with the fasta files prepared for phylogenetic analysis

conc=ProtNum + "_proteins_concatenation.fasta" # concatenation file name 

# declaring absolute file names
FullFile=basedir + conc
rate4site=r4sDir + "/Normalized_r4s_" + ProtNum + "_proteins.tsv"
shannon=shannonDir + "/Normalized_shannon_" + ProtNum + "_proteins.tsv"

data=FASTA_iterator(FullFile)

# coding protein names with letters
protdict={
	'A':'AHSG', 'B':'ALB', 'C':'AMBN', 'D':'AMELX', 'E':'AMTN', 'F':'COL17A1', 'G':'COL1A1',
	'H':'COL1A2', 'I':'COL2A1', 'J':'ENAM', 'K':'MMP20', 'L':'ODAM', 'M':'SERPINC1', 'N':'TUFT1'
}



def protein_counting(number):

	"""
	number accepted by this function is the index number of an amino acid position in a given protein concatenation.
	This function returns a letter, that represents the protein to which this amino acid belongs, code is defined in 
	dictionary "protdict"

	"""

	# protein coordinates relative to each concatenation, rg stands for range

	rg_14={
		1:[0,340], 2:[340,930], 3:[930,1347], 4:[1347,1525], 5:[1525,1716], 6:[1716,2701], 7:[2701,4142],
		8:[4142,5485], 9:[5485,6946], 10:[6946,8040], 11:[8040,8499], 12:[8499,8757], 13:[8757,9188], 14:[9188,9557],
	}

	rg_10={
		1:[0,340], 2:[340,930], 3:[930,1347], 4:[1347,1525], 5:[1525,1716], 6:[1716,2810], 
		7:[2810,3269], 8:[3269,3527], 9:[3527,3958], 10:[3958,4327], 
	}

	rg_5={
		1:[0,417], 2:[417,595], 3:[595,786], 4:[786,1880], 5:[1880,2339], 
	}


	if ProtNum == '14':
		if number in range(rg_14[1][0],rg_14[1][1]):
			return('A')
		elif number in range(rg_14[2][0],rg_14[2][1]):
			return('B')
		elif number in range(rg_14[3][0],rg_14[3][1]):
			return('C')
		elif number in range(rg_14[4][0],rg_14[4][1]):
			return('D')
		elif number in range(rg_14[5][0],rg_14[5][1]):
			return('E')
		elif number in range(rg_14[6][0],rg_14[6][1]):
			return('F')
		elif number in range(rg_14[7][0],rg_14[7][1]):
			return('G')
		elif number in range(rg_14[8][0],rg_14[8][1]):
			return('H')
		elif number in range(rg_14[9][0],rg_14[9][1]):
			return('I')
		elif number in range(rg_14[10][0],rg_14[10][1]):
			return('J')
		elif number in range(rg_14[11][0],rg_14[11][1]):
			return('K')
		elif number in range(rg_14[12][0],rg_14[12][1]):
			return('L')
		elif number in range(rg_14[13][0],rg_14[13][1]):
			return('M')
		elif number in range(rg_14[14][0],rg_14[14][1]):
			return('N')

	elif ProtNum == '10':

		if number in range(rg_10[1][0],rg_10[1][1]):
			return('A')
		elif number in range(rg_10[2][0],rg_10[2][1]):
			return('B')
		elif number in range(rg_10[3][0],rg_10[3][1]):
			return('C')
		elif number in range(rg_10[4][0],rg_10[4][1]):
			return('D')
		elif number in range(rg_10[5][0],rg_10[5][1]):
			return('E')
		elif number in range(rg_10[6][0],rg_10[6][1]):
			return('J')
		elif number in range(rg_10[7][0],rg_10[7][1]):
			return('K')
		elif number in range(rg_10[8][0],rg_10[8][1]):
			return('L')
		elif number in range(rg_10[9][0],rg_10[9][1]):
			return('M')
		elif number in range(rg_10[10][0],rg_10[10][1]):
			return('N')

	elif ProtNum == '5':
		
		if number in range(rg_5[1][0],rg_5[1][1]):
			return('C')
		elif number in range(rg_5[2][0],rg_5[2][1]):
			return('D')
		elif number in range(rg_5[3][0],rg_5[3][1]):
			return('E')
		elif number in range(rg_5[4][0],rg_5[4][1]):
			return('J')
		elif number in range(rg_5[5][0],rg_5[5][1]):
			return('K')

###########################################################################

# For Rate4site

print("Building rate4site subset fasta files")

fasta_outdir=basedir + "/subsets"

r4s_Un_file=fasta_outdir + "/rate4site/Under_" + ProtNum + "_proteins.fasta"
r4s_Ov_file=fasta_outdir + "/rate4site/Over_" + ProtNum + "_proteins.fasta"

for seq in data: # for each sequence in full concatenation fasta file
	
	# Under and Over threshold 1 of shannon and r4s
	under_list=[]
	over_list=[]
	
	# Lists with letter coding for the aa's belonging to each protein
	under_prot_list=[]
	over_prot_list=[]

	tot_len=len(seq[1])

	header=seq[0]
	with open(rate4site,"r") as r4s: # file with rate for site scorings per positon
		for line in r4s:

			line=line.strip()
			index=int(line.split("\t")[0])-1
			score=float(line.split("\t")[1])


			if score <= 1:

				under_prot_list.append(protein_counting(index)) # appending the code of protein to which this site belongs	
				under_list.append(seq[1][index]) # appending the aminoacid to the subset list for building the subset file
				
			else:

				over_prot_list.append(protein_counting(index)) # appending the code of protein to which this site belongs
				over_list.append(seq[1][index]) # appending the aminoacid to the subset list for building the subset file


	# Building the subset fasta files
	under_seq="".join(under_list)

	with open(r4s_Un_file, "a") as r4s_un:

		r4s_un.write(header+"\n"+ under_seq + "\n")

	over_seq="".join(over_list)

	with open(r4s_Ov_file, "a") as r4s_ov:

		r4s_ov.write(header+"\n"+ over_seq + "\n")
	

# Sequences of letter coding by protein for Rate4site
under_prot_seq="".join(under_prot_list)
over_prot_seq="".join(over_prot_list)

proteins=list(set(under_prot_seq))
proteins.sort()

###########################################################################

# For shannon entropy

print("Building shannon subset fasta files")

sh_Un_file=fasta_outdir + "/shannon/Under_" + ProtNum + "_proteins.fasta"
sh_Ov_file=fasta_outdir + "/shannon/Over_" + ProtNum + "_proteins.fasta"

data=FASTA_iterator(FullFile)

for seq in data: # for each sequence in full concatenation fasta file
	
	# Under and Over threshold 1 of shannon
	under_sh_list=[]
	over_sh_list=[]
	
	# Lists with letter coding for the aa's belonging to each protein
	under_prot_sh_list=[]
	over_prot_sh_list=[]

	tot_sh_len=len(seq[1])

	header=seq[0]
	with open(shannon,"r") as sh: # file with shannon entropy scorings per positon
		for sh_line in sh:

			sh_line=sh_line.strip()
			sh_index=int(sh_line.split("\t")[0])-1
			sh_score=float(sh_line.split("\t")[1])


			if sh_score <= 1:

				under_prot_sh_list.append(protein_counting(sh_index)) # appending the code of protein to which this site belongs	
				under_sh_list.append(seq[1][sh_index]) # appending the aminoacid to the subset list for building the subset file
				
			else:
				# print(sh_score)
				over_prot_sh_list.append(protein_counting(sh_index)) # appending the code of protein to which this site belongs
				over_sh_list.append(seq[1][sh_index]) # appending the aminoacid to the subset list for building the subset file


        # Building the subset fasta files
	under_sh_seq="".join(under_sh_list)

	with open(sh_Un_file, "a") as sh_un:

		sh_un.write(header+"\n"+under_sh_seq + "\n")

	over_sh_seq="".join(over_sh_list)

	with open(sh_Ov_file,"a") as sh_ov:

		sh_ov.write(header+"\n"+over_sh_seq + "\n")


# Sequences of letter coding by protein for shannon entropy
under_sh_prot_seq="".join(under_prot_sh_list)
over_sh_prot_seq="".join(over_prot_sh_list)

sh_proteins=list(set(under_sh_prot_seq))
sh_proteins.sort()

###########################################################################

# Building dictionaries for Rate4site and Shannon entropy coordinates both for Under, and Over values
# These dictionaries are used for building the IQtree nexus files

counter=0
UnResult=0
OvResult=0

UnResult_sh=0
OvResult_sh=0

r4sUnDict={}
r4sOvDict={}

shUnDict={}
shOvDict={}

for thyprot in proteins:
	UnCount=under_prot_seq.count(thyprot)
	OvCount=over_prot_seq.count(thyprot)

	UnCount_sh=under_sh_prot_seq.count(thyprot)
	OvCount_sh=over_sh_prot_seq.count(thyprot)

	if counter == 0:
		
		print(protdict[thyprot] + " = " + str(1) + "-" + str(UnCount))
		print(protdict[thyprot] + " = " + str(1) + "-" + str(OvCount))
		
		# range appended to dictionary is 1 to count of aminoacids over/under threshold of first protein in loop
		# for Rate4site
		r4sUnDict[protdict[thyprot]]= str(1)+"-"+str(UnCount)	
		r4sOvDict[protdict[thyprot]]= str(1)+"-"+str(OvCount)
		
		# for Shannon entropy
		shUnDict[protdict[thyprot]]= str(1)+"-"+str(UnCount_sh)	
		shOvDict[protdict[thyprot]]= str(1)+"-"+str(OvCount_sh)	

		counter+=1

		UnResult=UnResult+UnCount
		OvResult=OvResult+OvCount

		UnResult_sh=UnResult_sh+UnCount_sh
		OvResult_sh=OvResult_sh+OvCount_sh

	else:

		theUnSum=UnCount+UnResult
		theOvSum=OvCount+OvResult

		theUnSum_sh=UnCount_sh+UnResult_sh
		theOvSum_sh=OvCount_sh+OvResult_sh

		print(protdict[thyprot] + " = " + str(UnResult+1) + "-" + str(theUnSum))
		print(protdict[thyprot] + " = " + str(OvResult+1) + "-" + str(theOvSum))
		
		# range appended to dictionary is previous count + 1, to count of aminoacids over/under threshold of next protein in loop
		r4sUnDict[protdict[thyprot]]= str(UnResult+1)+"-"+str(theUnSum)
		r4sOvDict[protdict[thyprot]]= str(OvResult+1)+"-"+str(theOvSum)

		shUnDict[protdict[thyprot]]= str(UnResult_sh+1)+"-"+str(theUnSum_sh)
		shOvDict[protdict[thyprot]]= str(OvResult_sh+1)+"-"+str(theOvSum_sh)

		UnResult=UnCount+UnResult
		OvResult=OvResult+OvCount

		UnResult_sh=UnCount_sh+UnResult_sh
		OvResult_sh=OvCount_sh+OvResult_sh

########################################################################################

# Building the coordinate nexus files for IQtree

r4s_UnCoord_file= basedir + "/r4s_"+ProtNum+"_proteins_under.nexus"
r4s_OvCoord_file= basedir + "/r4s_"+ProtNum+"_proteins_over.nexus"
sh_UnCoord_file= basedir + "/sh_"+ProtNum+"_proteins_under.nexus"
sh_OvCoord_file= basedir + "/sh_"+ProtNum+"_proteins_over.nexus"


if ProtNum == "5":
	model="\n\tcharpartition 5_proteins = JTT+I+G4:AMBN, JTTDCMut+F+G4:AMELX, JTT+R2:AMTN, JTT+R3:ENAM, JTT+I+G4:MMP20;\n\nend;"
elif ProtNum == "10":
	model="\n\tcharpartition 10_proteins = JTT+R3:AHSG, JTTDCMut+G4:ALB, JTT+I+G4:AMBN, JTTDCMut+F+G4:AMELX, JTT+R2:AMTN, JTT+R3:ENAM, JTT+I+G4:MMP20, JTT+I+G4:ODAM, JTT+I+G4:SERPINC1, JTT+G4:TUFT1;\n\nend; "
elif ProtNum == "14":
	model="\n\tcharpartition 14_proteins = JTT+R3:AHSG, JTTDCMut+G4:ALB, JTT+I+G4:AMBN, JTTDCMut+F+G4:AMELX, JTT+R2:AMTN, JTT+R3:COL17A1, JTTDCMut+F+R3:COL1A1, JTTDCMut+F+R3:COL1A2, JTTDCMut+F+R3:COL2A1, JTT+R3:ENAM, JTT+I+G4:MMP20, JTT+I+G4:ODAM, JTT+I+G4:SERPINC1, JTT+G4:TUFT1;\n\nend;"


# Building nexus file for r4site Under threshold


with open(r4s_UnCoord_file, "w" ) as r4s_U_f:
	r4s_U_f.write("#nexus\n\nbegin sets;\n")

	for prot,values in r4sUnDict.items():
		r4s_U_f.write("\tcharset "+ prot + "=" + values + ";\n")

	r4s_U_f.write(model)


# Building nexus file for r4site Over threshold

with open(r4s_OvCoord_file, "w" ) as r4s_O_f:
	r4s_O_f.write("#nexus\n\nbegin sets;\n")

	for prot,values in r4sOvDict.items():
		r4s_O_f.write("\tcharset "+ prot + "=" + values + ";\n")

	r4s_O_f.write(model)


# Building nexus file for shannon Under threshold

with open(sh_UnCoord_file, "w" ) as sh_U_f:
	sh_U_f.write("#nexus\n\nbegin sets;\n")

	for k,v in shUnDict.items():
		sh_U_f.write("\tcharset "+k + "=" + v + ";\n")

	sh_U_f.write(model)


# Building nexus file for shannon Over threshold

with open(sh_OvCoord_file, "w")as sh_O_f:
	sh_O_f.write("#nexus\n\nbegin sets;\n")

	for k,v in shOvDict.items():
		sh_O_f.write("\tcharset "+k + "=" + v + ";\n")

	sh_O_f.write(model)
