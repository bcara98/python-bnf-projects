#!/bin/bash

#Description: This program uses 2 fastq files, a reference genome, and a Mills file as an input to generate a gene variant calling file (vcf file).
#			  This program requires bwa, samtools, bcftools installed as well as java 8 and to have the GATK.jar file instaled in the $HOME/~bin directory 
#			  (with $HOME/~bin added to environment path)

reads1=""	#stores the first read file location
reads2=""	#stores the second read file location
ref=""		#stores the reference file location
output="" 	#stores the output filepath
millsFile=""	#stores the millsfile location
realign=0	#determines if a GATK realignment needs to be done
gunzip=0	#determines if output should be gunzipped
v=0		#determines if verbose mode is active
index=0		#determines if the final bam file needs to be indexed
answer="y"	#helps in whenever we need to overwrite an existing output file
help=0		#determines if help page should be opened


#uses the flags below to the needed argumens such as filepaths and other optional flags and stores them in the variables above
while getopts "a:b:r:o:f:ezvih" flag;
do
	case "$flag" in
		a)
			reads1=$OPTARG
			;;
		b)	
			reads2=$OPTARG
			;;
		r)
			ref=$OPTARG
			;;
		o)
			output=$OPTARG
			;;
		f)
			millsFile=$OPTARG
			;;
		e)
			realign=1
			;;
		z)
			gunzip=1
			;;
		v)
			v=1
			;;
		i)
			index=1
			;;
		h)	
			help=1
			;;
		*)
			;;
	esac

done

#determines if help page should be opened and if so it is displayed and the function exits
if [ $help == 1 ];then
	echo "This is the helping manual for the genevcf command"
	echo "Usage instructions:"
	echo -e "\t-a read1 \"Provide the first read file\""
	echo -e "\t-b read2 \"Provide the second read file\""
	echo -e "\t-r ref \"Provide the reference genome file\""
	echo -e "\t-f Mills_File \"Provide the mills file\""
	echo -e "\t-o output \"Specify the output file with the vcf extension\, if vcf file already exists user will be prompted whenever to overwrite it or not\""
	echo -e "\t-e \"Optional flag whenever user wants to perform realignment with GATK\""
	echo -e "\t-i \"Optional flag whenever user wants to index the final bam file\""
	echo -e "\t-z \"Optional flag whenever user wants to gunzip the final vcf file\""
	echo -e "\t-v \"Optional flag to execute the command in verbose mode\""
	echo -e "\t-h \"Open the helping manual\""
	exit 0
fi

#checks whenever the user entered all the needed files and if the file paths are valid and exist if not an error message is shown and function exits
if ! [ -f $reads1 ] || [ ${#reads1} == 0 ];then
	echo "Error: The first read file with name '$reads1' does not exist!"
	echo "Make sure you have enetered the correct file path!"
	exit 1
fi

if ! [ -f $reads2 ] || [ ${#reads2} == 0 ];then
	echo "Error: The second read file with name '$reads2' does not exist!"
	echo "Make sure you have enetered the correct file path!"
	exit 1
fi

if ! [ -f $ref ] || [ ${#ref} == 0 ];then
	echo "Error: The reference genome file with name '$ref' does not exist!"
	echo "Make sure you have enetered the correct file path!"
	exit 1
fi

if ! [ -f $millsFile ] || [ ${#millsFile} == 0 ];then
	echo "Error: The mills file with name '$millsFile' does not exist!"
	echo "Make sure you have enetered the correct file path!"
	exit 1
fi

#check whenever a specified output name is given if not an error message is shown and function exits
if [ ${#output} == 0 ];then
	echo "Error: You must specify an output for the VCF file"
	exit 1
fi

#checks if outoupt needs to be gunziped and helps determine if appropriate output file already exists
if [ $gunzip == 1 ];then
	#if needs to be gunziped it checks if output.vcf.gz exists and if so a warning message is dispalyed asking whenever the user wants to overwrite it
	if [ -f "${output}.vcf.gz" ] && [ ${#output} != 0 ];then
		echo "Warning: Output VCF file with the name '$output' already exists!"
		read -p "Do you want to overwrite this VCF file (y/n):" answer
	fi
else
	#if does not need to be gunzipped it checks if output.vcf exists and if so a warning message is dispalyed asking whenever the user wants to overwrite it
	if [ -f "${output}.vcf" ] && [ ${#output} != 0 ];then
		echo "Warning: Output VCF file with the name '$output' already exists!"
		read -p "Do you want to overwrite this VCF file (y/n):" answer
	fi
fi

#if the user does not want to overwrite the existing output file the function exits
if [ "$answer" != "y" ] && [ "$answer" != "Y" ];then
	exit 0
fi


#helper variables
dict_ref=$(echo "$ref"|sed 's/\(.*\)\(\..*$\)/\1.dict/') #determines the respecitive dictionary file name and path by removing the .fa extension and replaces it with .dict
cut_ref=$(echo "$ref"|sed 's/\(.*\/\)\(.*$\)/\2/') #only gets the reference genome filename without the full path and uses that to helpe generate a unique name for the intermediate bam and sam files


#For the reamining code below if verbose mode is on with v equal to 1 then it prints the steps the programm is executing. Bceause of this i made all the bwa, samtools, GATK, and bcftools supress their console output
if [ $v == 1 ];then
	echo "Indexing reference genome"
fi
#index the reference genome
bwa index $ref 2>/dev/null
if [ $v == 1 ];then
	echo "Complete"
	echo "Creating lane (this may take a while)"
fi
#generate the initial lane
bwa mem -R '@RG\tID:foo\tSM:bar\tLB:library1' $ref $reads1 $reads2 > ${cut_ref}_lane.sam 2>/dev/null
if [ $v == 1 ];then
	echo "Complete"
	echo "Fixing lane"
fi
#fix the lane with samtools
samtools fixmate -O bam ${cut_ref}_lane.sam ${cut_ref}_lane_fixmate.bam 2>/dev/null
if [ $v == 1 ];then
	echo "Complete"
	echo "Sorting lane"
fi
#sort the lane with samtools
samtools sort -O bam -o ${cut_ref}_lane_sorted.bam -T /tmp/lane_temp ${cut_ref}_lane_fixmate.bam 2>/dev/null
if [ $v == 1 ];then
	echo "Complete"
	echo "Generating fai index file for reference genome"
fi
#generate a fai index file for the reference genome with samtools
samtools faidx $ref 2>/dev/null
if [ $v == 1 ];then
	echo "Complete"
	echo "Generating dictionary file for reference genome"
fi
#generate a dictionary file for the reference genome
samtools dict $ref -o $dict_ref 2>/dev/null
if [ $v == 1 ];then
	echo "Complete"
fi
if [ $realign == 0 ];then
#if the user did not akd to realign then the program skips the GATK steps
	if [ $index == 1 ];then
	#if the user asked to index the final bam file which is just the sorted lane in this case else this step is skipped
		if [ $v == 1 ];then
			echo "Indexing Sorted lane"
		fi
		#sorted bam file is indexed with samtools
		samtools index ${cut_ref}_lane_sorted.bam 2>/dev/null	
		if [ $v == 1 ];then
			echo "Complete"
		fi
	fi
	if [ $gunzip == 1 ];then
	#if user used the -z flag to gunzip the vcf.gz extension is appended to the output filename
		if [ $v == 1 ];then
			echo "Generating gunziped VCF file"
		fi
		output="${output}.vcf.gz"
	else
	#if no gunzip is needed the vcf extension is appended to the output filename
		if [ $v == 1 ];then
			echo "Generating VCF file"
		fi
		output="${output}.vcf"
	fi
	#generate the vcf file gunzipped or not gunzipped depeniding on whenever the -z flag was used
	bcftools mpileup -Ou -f $ref ${cut_ref}_lane_sorted.bam 2>/dev/null| bcftools call -vmO z -o $output 2>/dev/null
	if [ $v == 1 ];then
		echo "Complete"
	fi
else
#if the user prompted to realign with the -e flag the the following statment is executed and realignment is performed with GATK
	if [ $v == 1 ];then
		echo "Indexing Sorted lane"
	fi
	#sorted lane is sorted regardless since it is requred by the GATK program
	samtools index ${cut_ref}_lane_sorted.bam 2>/dev/null	
	if [ $v == 1 ];then
		echo "Complete"
		echo "Generating the lane intervals"
	fi
	#the lane intervals are generated with the GATK program while it generates a logfile with the same name of the ouptut with the .log extension appeneded
	java -Xmx2g -jar $HOME/bin/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $ref -I ${cut_ref}_lane_sorted.bam -o ${cut_ref}_lane.intervals --known $millsFile &>"${output}.log"
	if [ $v == 1 ];then
		echo "Complete"
		echo "Realigning lane"
	fi
	#the sorted lane is realigned and any logs are appended to the previosuly creaded log file	
	java -Xmx4g -jar $HOME/bin/GenomeAnalysisTK.jar -T IndelRealigner -R $ref -I ${cut_ref}_lane_sorted.bam -targetIntervals ${cut_ref}_lane.intervals -known $millsFile -o ${cut_ref}_lane_realigned.bam &>>"${output}.log"
	if [ $v == 1 ];then
		echo "Complete"
	fi
	if [ $index == 1 ];then
	#if the indiex flag -i is used the the realigned bam file is indexed with samtools else it skips this step
		if [ $v == 1 ];then
			echo "Indexing realigned lane"
		fi
		samtools index ${cut_ref}_lane_realigned.bam 2>/dev/null
		if [ $v == 1 ];then
			echo "Complete"
		fi
	fi
	if [ $gunzip == 1 ];then
	#if user used the -z flag to gunzip the vcf.gz extension is appended to the output filename
		if [ $v == 1 ];then
			echo "Generating gunziped VCF file"
		fi
		output="${output}.vcf.gz"
	else
	#if no gunzip is needed the vcf extension is appended to the output filename
		if [ $v == 1 ];then
			echo "Generating VCF file"
		fi
		output="${output}.vcf"
	fi
	#generate the vcf file gunzipped or not gunzipped depeniding on whenever the -z flag was used
	bcftools mpileup -Ou -f $ref ${cut_ref}_lane_realigned.bam 2>/dev/null| bcftools call -vmO z -o $output 2>/dev/null
	if [ $v == 1 ];then
		echo "Complete"
	fi
fi

#deleting all the intermediate files
if [ $v == 1 ];then
	echo "Deleting intermediate files"
fi
rm -f ${cut_ref}*.bam
rm -f ${cut_ref}*.bai
rm -f ${cut_ref}*.sam
rm -f ${cut_ref}*.intervals
rm -f ${ref}.*
rm -f $dict_ref
rm -f ${millsFile}.idx

if [ $v == 1 ];then
	echo "All steps completed succesfully"
fi

