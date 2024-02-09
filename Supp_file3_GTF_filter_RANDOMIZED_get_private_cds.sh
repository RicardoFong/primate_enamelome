#!/bin/bash


###################################################################################################################
###################################################################################################################

# autor: Ricardo Fong Zazueta
# mantainer: Ricardo Fong Zazueta
# mail: ricardo.fong.zazueta@umontreal.ca

# This code outputs the two alelles when the gene is heterozygous, and just one alelle when it is homozygous. In the 
# latter case, the option used calls all the alternative nucleotide when a snp is found.

# The code not yet makes any filtering of the vcf variants, therefore the vcf data used should be previously filtered
# in order to get biologically relevant data. Working on implementing the filters.

# Usage:

# ./thisprogram.sh gtf_file_location.gtf ref_genome_location.fa vcf_file.vcf sample_name output_absolute_dir Species_name

# Three files are going to be written to the directory you give as the fifth argument in the "Usage". Two fasta files,
# one with the protein sequence of the input genes, and one with the correspondent DNA sequence. The third file is an intermediate
# tsv extension file with the CDS's retreived from the GTF file.

# IMPORTANT: The version of this program as of the 3th dec 2021 (RFZ_local_get_private_cds.sh) expects the name
# of the species whith the one you are working as a sixth argument in the initial command. This name will be appended
# to the fasta headers in the ouput.

# In order to translate proteins, the correspondent RFZ_local_combine_cds_to_prot.py file must be in the folder where
# you are running this bash program.


# As stated in the usage, the fifth is the output directory. It is very important for it to be an absolute path, or a path
# that is contained in your current working directory.

###################################################################################################################
###################################################################################################################

module load gcc/6.3.0
module load xz/5.2.2
module load BCFTOOLS/1.9
module load SAMTOOLS/1.9
module load zlib/1.2.4
module load PYTHON/3.8.6
module load TABIX/0.2.6

echo

gtf=$1
fasta=$2
mainvcf=$3
sample=$4
base_output=$5 
species=$6

gene_ids=$(cat $gtf | awk '{print $10}'  | sed -e 's/\;//g' -e 's/\"//g' | sort | uniq)


# Getting the sample name
# If a simple sample name is provided, it is respected as-it-is

rm -f ${base_output}/${id}.private_exon_sequences.tsv 2> /dev/null


for identifier in ${gene_ids[@]}; do

	# Checking if the gene in the loop is heterozygous or not

	Het_check=$( cat $gtf | grep $identifier | grep CDS | sort -nk4 | \
	while read l;do 
    
		scaffold=$(echo $l | awk '{print $1}');
		start_co=$(echo $l| awk '{print $4}');
		stop_co=$(echo $l | awk '{print $5}');

		snpgeno=$(bcftools view $mainvcf --samples $sample ${scaffold}:${start_co}-${stop_co} | \
		bcftools filter -e "TYPE!='snp' | (GT='het' & FMT/AD[*:*] < 3 ) | AC > 2 | FMT/DP <= 10 | \
		QD < 2 | FS > 60 | INFO/MQ < 40 | SOR > 3 | ReadPosRankSum < -8 | QUAL < 30 | MQRankSum < -12.5" | \
		grep -v "^#" | awk '{print $10}' | cut -d ':' -f 1)
			
		indelgeno=$(bcftools view $mainvcf --samples $sample ${scaffold}:${start_co}-${stop_co} | \
                bcftools filter -e "TYPE!='snp' | (GT='het' & FMT/AD[*:*] < 3 ) | FMT/DP <= 10 | \
                QD < 2 | FS > 200 | INFO/MQ < 40 | SOR > 3 | ReadPosRankSum < -20 | QUAL < 30" | \
                grep -v "^#" | awk '{print $10}' | cut -d ':' -f 1)
		

		for i in $snpgeno; do if [[ $i == "0/1" ]]; then echo 1;
			 
			fi; 
		done;
		
		for i in $indelgeno; do if [[ $i == "0/1" ]]; then echo 1;

			fi;
		done;
	done; )
		
	if [ -z "$Het_check" ]; then

#############################################################################################################

		# This code is run if gene is homozygous
		echo $identifier
		echo $sample is homozygous for $identifier, getting private exons > /dev/stderr
 		
		cat $gtf | grep $identifier | grep CDS | sort -nk4 | while read l;do 
    			
			scaffold=$(echo $l | awk '{print $1}');
			start_co=$(echo $l| awk '{print $4}');
			stop_co=$(echo $l | awk '{print $5}');
			orientation=$(echo $l | awk '{print $7}');
			gene=$(echo $l | awk '{print $10}' | sed  -e 's/\;//g' -e 's/\"//g');
			tag=$(echo $l|awk '{print $12}' | sed -e 's/\;//g' -e 's/\"//g')"_A1";
			
			bcftools view $mainvcf --samples $sample ${scaffold}:${start_co}-${stop_co} | \
	                bcftools filter -e "TYPE!='snp' | (GT='het' & FMT/AD[*:*] < 3 ) | AC > 2 | FMT/DP <= 10 | \
        	        QD < 2 | FS > 60 | INFO/MQ < 40 | SOR > 3 | ReadPosRankSum < -8 | QUAL < 30 | MQRankSum < -12.5" > "tmp.snp.vcf"
			
			bgzip -f tmp.snp.vcf
			tabix -f tmp.snp.vcf.gz

			bcftools view $mainvcf --samples $sample ${scaffold}:${start_co}-${stop_co} | \
                        bcftools filter -e "TYPE!='indel' | (GT='het' & FMT/AD[*:*] < 3 ) | FMT/DP <= 10 | \
                        QD < 2 | FS > 200 | INFO/MQ < 40 | SOR > 3 | ReadPosRankSum < -20 | QUAL < 30" > "tmp.indel.vcf"
			
			bgzip -f tmp.indel.vcf
			tabix -f tmp.indel.vcf.gz

			samtools faidx $fasta ${scaffold}:${start_co}-${stop_co} | \
			bcftools consensus --sample $sample -i "GT='AA'" tmp.snp.vcf.gz | \
			bcftools consensus --sample $sample -i "GT='AA'" tmp.indel.vcf.gz | \
			/scratch/devel/jkrueger/bin/bioawk-master/bioawk -cfastx \
				-vtag=$tag \
				-vstop_co=$stop_co \
				-vstart_co=$start_co \
				-vscaffold=$scaffold\
				-vorientation=$orientation\
				-vgene=$gene\
			'{OFS="\t";print scaffold, start_co-1, stop_co, orientation, tag, $seq, gene}' | \
			sort -k1,1 -k2,2n >> ${base_output}/$id.private_exon_sequences.tsv;
		done;

	else

#############################################################################################################

		# This code is run if gene is heterozygous
		rand=$(perl therandom.pl)
		
		if [ $rand -ge 4 ]; then exp=1; else exp=2; fi; 
		
		echo $sample is heterozygous for $identifier, getting private exons > /dev/stderr
  
		cat $gtf | grep $identifier | grep CDS | sort -nk4 | while read l;do 
    
			scaffold=$(echo $l | awk '{print $1}');
			start_co=$(echo $l| awk '{print $4}');
			stop_co=$(echo $l | awk '{print $5}');
			orientation=$(echo $l | awk '{print $7}');
			gene=$(echo $l | awk '{print $10}' | sed  -e 's/\;//g' -e 's/\"//g');
			tag=$(echo $l|awk '{print $12}' | sed -e 's/\;//g' -e 's/\"//g')"_A"$exp;

			bcftools view $mainvcf --samples $sample ${scaffold}:${start_co}-${stop_co} | \
	               	bcftools filter -e "TYPE!='snp' | (GT='het' & FMT/AD[*:*] < 3 ) | AC > 2 | FMT/DP <= 10 | \
        	       	QD < 2 | FS > 60 | INFO/MQ < 40 | SOR > 3 | ReadPosRankSum < -8 | QUAL < 30 | MQRankSum < -12.5" > "tmp.snp.vcf"

			bgzip -f tmp.snp.vcf
			tabix -f tmp.snp.vcf.gz
			
			bcftools view $mainvcf --samples $sample ${scaffold}:${start_co}-${stop_co} | \
                       	bcftools filter -e "TYPE!='indel' | (GT='het' & FMT/AD[*:*] < 3 ) | FMT/DP <= 10 | \
                       	QD < 2 | FS > 200 | INFO/MQ < 40 | SOR > 3 | ReadPosRankSum < -20 | QUAL < 30" > "tmp.indel.vcf"
			
			bgzip -f tmp.indel.vcf
			tabix -f tmp.indel.vcf.gz

			samtools faidx $fasta ${scaffold}:${start_co}-${stop_co} | \
			bcftools consensus --sample $sample -H $exp tmp.snp.vcf.gz | \
			bcftools consensus --sample $sample -H $exp tmp.indel.vcf.gz | \
			/scratch/devel/jkrueger/bin/bioawk-master/bioawk -cfastx \
				-vtag=$tag \
				-vstop_co=$stop_co \
				-vstart_co=$start_co \
				-vscaffold=$scaffold\
				-vorientation=$orientation\
				-vgene=$gene\
			'{OFS="\t";print scaffold, start_co-1, stop_co, orientation, tag, $seq, gene}' | \
			sort -k1,1 -k2,2n >> ${base_output}/$id.private_exon_sequences.tsv;
		done;
		
#############################################################################################################

	fi;

done;

echo

echo combining exons of provided gtf, genes translated are: $gene_ids > /dev/stderr

python3 Combine_CDS_concatenator.py ${base_output}/$id.private_exon_sequences.tsv $id $species

rm tmp.snp.vcf.gz* tmp.indel.vcf.gz*
