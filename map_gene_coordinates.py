#!/usr/bin/env python3

#Description: This program maps the gene ID with its coordinates from the UCSC genome browser to the gene name in the kgXref.txt file and from that it filters it
#             based on a third file (like list of infectious gene names)into an output file the gene name, its corresponding chromosome and coordinates 

import sys

#extracts file paths from command line arguments
known_gene_path = sys.argv[1]
kgxref_path = sys.argv[2]
infectious_sets_path = sys.argv[3]
#opens each file in read mode stores their line contents into lists and then closes the files
known_gene = open(known_gene_path,'r')
known_gene_content = known_gene.readlines()
known_gene.close()
kgxref = open(kgxref_path,'r')
kgxref_content = kgxref.readlines()
kgxref.close()
infectious_sets = open(infectious_sets_path,'r')
infectious_sets_content = infectious_sets.readlines()
infectious_sets.close()
#delete the first line read from the infectious sets file that only has the header "Infectious Diseases"
del infectious_sets_content[0]

#dictionaries uses to speed up indexing process into constant O(1) time and to easier prevent duplicates
genes = {} #uses the gene name from the kgxref file as the key and the unique UCSC gene ID as the value (only saves first instance of the gene found, ignored other instances)
gene_id = {} #uses the unique UCSC gene ID as the key from the knownGene file and a list with the chromosome number start and stop positions as the value

for i in range (0,len(known_gene_content)):
    #remove newline characters from the list containing the knownGene file contents
    known_gene_content[i] = known_gene_content[i].replace('\n',"")
    known_gene_content[i] = known_gene_content[i].split('\t') #splits each line based on the tab character
    #adds element into the gene_id dictionary (key: UCSC gene ID, value: list of chromosome number start and stop positions)
    gene_id[known_gene_content[i][0]] = [known_gene_content[i][1],known_gene_content[i][3],known_gene_content[i][4]]
for i in range (0,len(kgxref_content)):
    #remove newline characters from the list containing the kgxref file contents
    kgxref_content[i] = kgxref_content[i].replace('\n',"")
    kgxref_content[i] = kgxref_content[i].split('\t') #splits each line based on the tab character
    #if gene dictionary already has a gene name skip this step else add to the dictionary with this format (key: gene name, value: UCSC gene ID)
    if kgxref_content[i][4] in genes:
        pass
    else:
        genes[kgxref_content[i][4]] = kgxref_content[i][0]
for i in range (0,len(infectious_sets_content)):
    #remove newline characters from the list containing the infectious sets file contents
    infectious_sets_content[i] = infectious_sets_content[i].replace('\n',"")
    infectious_sets_content[i] = infectious_sets_content[i].split('\t')[0] #gets the infectious gene name from each line


infectious_sets_content.sort() #sort the infectious gene set names in an alphabetically sorted manner
output = open("geneSetCoordinates.txt",'w') #create new file to save output contents
output.write('Gene\tChr\tStart\tStop\n') #write header into file
print('Gene\tChr\tStart\tStop') #print header into console
for i in infectious_sets_content:
    #if gene name in infectious set file exists in kgxref file print output and also add it to the output file else skip gene
    if i in genes:
        output.write(i+'\t'+gene_id[genes[i]][0]+'\t'+gene_id[genes[i]][1]+'\t'+gene_id[genes[i]][2]+'\n')
        print(i,'\t',gene_id[genes[i]][0],'\t',gene_id[genes[i]][1],'\t',gene_id[genes[i]][2])
output.close() #closes ouput file
