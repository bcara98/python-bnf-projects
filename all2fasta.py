#!/usr/bin/env python3 
import argparse
import re

#Description: This program converts any relevant sequence from a varienty of file types to a fasta format. The extension of the input files does not matter as the program
#             scans for content instead. The program also allows you to enter the fold number.

#parses the fold and input file arguments
parser = argparse.ArgumentParser(description="Convert any sequence relevant file type to fasta format")
parser.add_argument('-f',metavar='--fold',type=int,default=70,help='Enter the fold number (How many nucleotides/amino acids to display per line-Default 70)')
parser.add_argument('-i',metavar='--input',type=str,help='Enter the input file to convert to fasta format')

args = parser.parse_args()

#helper function to check if either arguments are valid by passing a list of boolean values
def arg_check(valid):
    for i in valid:
        if i == False:
            print("Terminating program!") #if one of the arguments is invalid ther the program terminates
            exit(1)

#helperf function to detect the file type by scanning the file contents with the help of regex expressions
def detect_file_type(content):
    #regex expression for fastq content including the @ followed by any non white space character
    fastq_pattern = r'^@[^\s]*\n'
    #regex expression for sam content that includes a bunch of lines starting with @ and at least one of them is followed by the SQ and PG flag 
    sam_pattern = r'^(@[^\s].*\n)*(@SQ.*\n)(@[^\s].*\n)*(@PG.*\n)|^(@[^\s].*\n)*(@PG.*\n)(@[^\s].*\n)*(@SQ.*\n)'
    #regex expression for mega content scanning for the mega header
    mega_pattern = r'^#MEGA\s*|^#mega\s*\n'
    #regex expression for embl content scanning for the ID label at the begining of the file
    embl_pattern = r'^ID\s*.*\n'
    #regex expression for genebank content scanning for the LOCUS label at the begining of the file
    genebank_pattern = r'^LOCUS\s*.*\n'
    #regex expression for vcf content scanning for the ##fileformat=VCF header at the begining of the file
    vcf_pattern = r'^##fileformat=VCF.*\n'
    #store regex patterns in a list
    patterns = [fastq_pattern,sam_pattern,mega_pattern,embl_pattern,genebank_pattern,vcf_pattern]
    j = -1 
    #loop scans for match in patterns list, if found j will have the index value of that last with the matched pattern
    for i in range(0,len(patterns)):
        if re.match(patterns[i],content):
            j = i

    #switch statement based on the value of j derived from the patterns list
    switch={
        0:'fastq',
        1:'sam',
        2:'mega',
        3:'embl',
        4:'genebank',
        5:'vcf'
    }
    file_type = switch.get(j) #this detects the appropriate file type
    return file_type #returns the file type as a string

#helper function that helps scan though the file contonet and extract the sequences along with their names based on the file type
def parse_content(content,file_type):
    seq = {} #dictionary stores sequence IDs as key and their corresponding sequence as value

    #if the file type is fastq
    if file_type == 'fastq':
        j = 0 #keeps track of the mutiple of each line number since you can have 
        seq_id = ''
        for i in content:
            if len(i) != 0: #if string is not empty then check following
                if j == 0: 
                    seq_id = i[1:] #if it is the first line of the entry (sequence name) just get the sequence name and store it as key to the dictionary
                if j == 1: 
                    seq[seq_id] = i #if it is the second line in the entry (the sequence) just add it to the dictionary as a vlaue with its corresponding name as key
                j += 1 #incrment line multiple
                if j == 4: 
                    j = 0 #if reaching the end for one entry reset j to 0
    #if the file type is embl
    elif file_type == 'embl':
        seq_id = '' #stores the sequence name
        sequence = '' #stores the sequence ID
        found_sequence = False #will become true if sequence found

        #loop scans thrhough file
        for i in content:
            if re.match(r'^ID\s*.*',i):
                #when are at the ID section it extracts the sequence name with the version with the help of regex expressions and process them accordingly to genereate the name in the proper format
                seq_id += i 
                seq_id = re.sub(r'([^\s]*\s*)([A-Za-z0-9]*)(;\s*)(SV\s*)([0-9]*)(;.*)',r'ENA|\2|\2.\5 ',seq_id)
            if re.match(r'^DE\s*.*',i):
                seq_id += re.sub(r'(DE\s*)(.*)',r'\2',i) #when at DE section additionaly info is added to the sequence ID
            if re.match(r'^SQ\s*',i): 
                found_sequence = True #if reached at SQ section you are at sequence so this value becomes true
                continue
            if found_sequence == True: 
                sequence+=i #if sequence is found then it is added to the sequence string

        sequence = re.sub(r'[^A-Za-z]','',sequence) #regex substitution helps to remove any non alphabetic characters (numerical, white space)
        seq[seq_id] = sequence #adds the sequence ID and seequence to dictionary
    #if the file type is genebank
    elif file_type == 'genebank':
        seq_id = '' #stores the sequence name
        sequence = '' #stores the sequence ID
        found_sequence = False #will become true if sequence found

        #loop scans thrhough file
        for i in content:
            if re.match(r'^VERSION\s*.*',i):
                #when are at the VERSION section it extracts the sequence version with help of regex expressions and adds it at the begining of the seq_id string to genereate the name in the proper format
                seq_id = re.sub(r'(VERSION\s*)([^\s]*)',r'\2',i) + ' ' + seq_id
            if re.match(r'^DEFINITION\s*.*',i):
                #when are at the DEFINTION section it extracts additional info about the sequence name
                seq_id += re.sub(r'^(DEFINITION\s*)(.*)',r'\2',i)
            if re.match(r'^ORIGIN\s*',i):
                #if reached at ORIGIN section the you are at the sequence so value becomes true
                found_sequence = True
                continue
            if found_sequence == True:
                sequence+=i #if sequence is found then it is added to the sequence string

        sequence = re.sub(r'[^A-Za-z]','',sequence) #regex substitution helps to remove any non alphabetic characters (numerical, white space)
        seq[seq_id] = sequence #adds the sequence ID and seequence to dictionary
    #if file type is mega
    elif file_type == 'mega':
        passed_header = False #checks if you passed mega header
        found_title = False #checks if found sequence name for sequential mega format
        seq_id = ''
        sequence = ''
        for i in content:
            if re.match(r'^#MEGA\s*|^#mega\s*',i) and passed_header == False:
                #detects if you passed mega header and sets value to true
                passed_header = True
                continue
            if re.match(r'^#[^\s]*\s+[^\s]+',i) and passed_header == True:
                #scans for sequences if mega file is in an interleaved format
                seq_id = re.sub(r'^(#)([^\s]*)(\s*)([^\s]*)(\s*)', r'\2',i) #gets sequence ID though regex
                if seq_id in seq:
                    seq[seq_id] += re.sub(r'^(#)([^\s]*)(\s*)([^\s]*)(\s*)', r'\4',i) #appends sequence portion to seq dictionary if extracted seq_id already exists
                else:
                    seq[seq_id] = re.sub(r'^(#)([^\s]*)(\s*)([^\s]*)(\s*)', r'\4',i) #if seq_id is not in dictionary it adds it with the first found sequence portion
                continue
            if re.match(r'^#[^\s]*',i) and passed_header == True:
                #detects if you found the title (always after passing the mega header) for sequeential mega format
                seq_id = re.sub(r'^(#)([^\s]*)', r'\2',i) #adds title to seq_id string
                found_title = True
                continue
            if found_title == True:
                #scans sequence and adds it to dictionary in sequential mega format
                if seq_id in seq:
                    seq[seq_id] += i
                else:
                    seq[seq_id] = i
            if i == '':
                found_title = False #reset found title to false when reaching empyt line (about to go to next sequence in sequential mega format)
    #if file type is sam
    elif file_type == 'sam':
        #scans through sam file and gets the sequence ID from the first column and sequence from the 10th column with help of regex and addes them to dictionary accordingly
        for i in content:
            if re.match(r'^[^@].*',i):
                seq_id = re.sub(r'^([^\s]*)(\s*.*)',r'\1',i)
                sequence = re.sub(r'^([^\s]*\s*){9}([^\s]*)(\s*.*)',r'\2',i)
                seq[seq_id] = sequence
    #if file type is vcf
    elif file_type == 'vcf':
        sample_list = [] #stores the number of samples (including extra entry for reference genome)
        ref_sequence = [] #stores the whole reference sequence in segments
        alt_sequence = [] #stores the whole alternate sequence in segments
        sample_info = [] #stores info about samples (which portions refere to reference and alternate sequence. It is a 2D list)
        added_seq = False #checks if extra entry for reference sample has been added to sample list later in the loop

        #loop scans the contents
        for i in content:
            if re.match(r'^#CHROM.*',i):
                #if you are at column titles it extracts all samples and separates them based on the tab-delimeter
                i = re.sub(r'(^#CHROM.*FORMAT\t)(.*)',r'\2',i)
                sample_list = i.split('\t')
                continue
            if re.match(r'^[^#].*',i):
                #first entry after the column titles should contain the reference sample name and sets added_seq to true
                added_seq = True
                seq_id = re.sub(r'^([^\s]*)(\s*.*)',r'\1',i) #regex extacts the sequence name
                if not (added_seq == True and seq_id == sample_list[0]):
                    sample_list.insert(0,seq_id) #if reference sequence name is not added to sample list then add it at begining of list

                #scans vertically in the file through each iteration to get the reference and alternate sequence segemtns and store them to their respecive lists
                ref_sequence.append(re.sub(r'^([^\s]*\s*){3}([^\s]*)(\s*.*)',r'\2',i))
                alt_sequence.append(re.sub(r'^([^\s]*\s*){4}([^\s]*)(\s*.*)',r'\2',i))
                l = [] #will get each line that stores info about each sample at that specific location if they use the reference or alternate sequence segment
                for j in range(1,len(sample_list)):
                    re_str = '^([^\s]*\s*){'+str(j+8)+'}([0-9]*)(:\s*.*)' #regex string will be modified at every looop iteration to scan horizontialy amoing all the samples
                    l.append(re.sub(re_str,r'\2',i)) #the needed info is added to temporary list l
                sample_info.append(l) #list l is then appended to the 2D sample list
        #there could be multiple alternative sequences if so then split them based on the ',' character
        for i in range(0,len(alt_sequence)):
            alt_sequence[i] = alt_sequence[i].split(',')
        #extract the reference sequence and pair it with the reference sample name as you add it to dictionary accordingly
        ref_seq = '' 
        for i in ref_sequence:
            ref_seq+=i
        seq[sample_list[0]] = ref_seq
        #nested for lop scans verically first for each sample and then moves horizationaly
        for i in range(0,len(sample_info[0])):
            seq_str = '' #set sequence to empty string
            for j in range(0,len(sample_info)):
                #scan vertically for each sample
                if int(sample_info[j][i]) == 0:
                    seq_str += ref_sequence[j] #if value is 0 in sample list then refer to reference sequence for that particular segment and add it to seq_str
                else:
                    seq_str += alt_sequence[j][int(sample_info[j][i])-1] #if value is 0 in sample list then refer to appropriate alternate sequence for that particular segment and add it to seq_str
            seq[sample_list[i+1]] = seq_str #add the sample name and correctly extracted sequence to dictionary
    return seq #returns the dictionary with the sequence names and respective sequences

#scans the sequences in the dictionary to detect if it is an amino acid or nucleotided sequence
def get_seq_type(seq_dict):
    for i in seq_dict.values():
        if re.match(r'^.*[^ACGTNacgtn]+.*$',i): 
            return 'P' #return 'P' if any other character besides ACGTNactn is found in the sequence
    return 'N' #else return 'N'

#helper function to generate the fasta file with the correct format
def generate_fasta_string(seq_dict,seq_type):

    #generates the fasta file name, if sequence type is amino acid it adds the .faa extension else it addes the .fna extension
    file_name = re.sub(r'(.*[^\.]+)(\..*)$',r'\1',args.i)
    if seq_type == 'P':
        file_name+='.faa'
    else:
        file_name+='.fna'

    out = open(file_name,'w') #creates output file based on the correctly generated file name
    #loop scans for all sequences in the dictionary        
    for i in seq_dict.keys():
        head = '>'+i+'\n' #string will store the header file with the sequence name
        seq = '' #will store respective sequence correctly folded
        k = 0 #keeps track of the fold

        #scans each character in sequence
        for j in seq_dict[i]:
            #if character count is larger than the provided fold number then add a new line and reset k to 0
            if k >= args.f:
                k = 0
                seq += '\n'
            seq+=j
            k+=1
        seq = seq.upper() #make sure every character in  sequence is uppercase
        out.write(head+seq+'\n') #write the header with the sequence to the file

valid = [True,True] #keeps track if arguments are valid
if args.f < 1:
    #if fold number is invalid
    print("Invalid fold number, must be larger than 0!")
    valid[0] = False

#tries to open file if not able to open (does not exists it throws exception)
try:
    file = open(args.i,'r')
except IOError as x:
    print("Could not open file, or file does not exist!")
    valid[1] = False

arg_check(valid) #passes the valid list if one of them is false terminates program
      
content = file.read() #reads file contents stores them to content
file.close() #closes input file
file_type = detect_file_type(content) #detects file type based on content
content = content.split('\n') #split content into list based on the new line characters
seq_list = parse_content(content,file_type) #properly parse content as a list and file type to extract the sequences and their names
t = get_seq_type(seq_list) #gest the sequence type (amino acid, nucleotide)
generate_fasta_string(seq_list,t) #generates the fasta file based on sequence type and dictionary with sequences and their names
