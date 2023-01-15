#!/usr/bin/env python3 
import argparse
import re

#Description: This program counts the number of occurences of each chromosme region in a BED file. 
#             Ex: [chr1 10 15] and [chr1 12 17], the follwing regions result (chr1 10,12->1 time, chr1 12 15->2 times, chr1 15 17->1 time)

#parses the input file argument
parser = argparse.ArgumentParser(description="Counts the occurences each chromosome region is referenced in the BED file")
parser.add_argument('-i',metavar='--input',type=str,help='Enter the input BED file')

args = parser.parse_args()

#helper function to help determine the number of occurences in each chromosome region
def parse_occurences(content):
    chr_pos = {} #dictionary will stores each chromomsome along with a list of checkpoints to scan through
    #loop scans through the file content
    for i in content:
        #if chromosome is not in dictionary add it as key and an empty list as value
        if i[0] not in chr_pos:
            chr_pos[i[0]] = []
        #appends the starting and ending points it reades for each chromosome to respective list in dictionary
        #from bed file as a list of 2 elemens [Position number,S=Start/E=End label]
        chr_pos[i[0]].append([int(i[1]),'S'])
        chr_pos[i[0]].append([int(i[2]),'E'])
    
    #scans each chromosome key in dictionary and sorts elements of the list it conatins based on the first value of the sublist (Positon number)
    for i in chr_pos.keys():
        chr_pos[i].sort()

    occurence_list = [] #stores sublists in this format [chromosome Number, start postion, end position, coverage] it is a 2D list
    occurence = 0
    for i in chr_pos.keys():
        l = [] #temporary 2D helper list to scan for each region with a specific occurence for each chromosome, resets when scanning different chromosome
        for j in chr_pos[i]:
            #while scanning each crhomomosome

            if j[1] == 'S':
                occurence += 1 #if position number is a starting position increment occurences by 1
            elif j[1] == 'E':
                occurence -= 1 #if position number is an ending position decrement occurences by 1
            l.append([j[0],occurence]) #add the corresponding position number and occurence to temporary list l
            #if list l stores 2 points
            if len(l) >= 2:
                #add they are neither the same postion and neither occurence in first point is 0
                if l[0][0] != l[1][0] and l[0][1] != 0:
                    #add to the occurence list a list with chromosome number, start position, end postion, coverage (occurence) number for that region
                    occurence_list.append([i,l[0][0],l[1][0],l[0][1]]) 
                l = [l[1]] #reset l to only hold the second element
    return occurence_list #returns the list containing the coverage info for all the chromosomes in the bed file 

#tries to open file if not able to open (does not exists it throws exception)           
try:
    file = open(args.i,'r')
except IOError as x:
    print("Could not open file, or file does not exist!")
    exit(1)

content = file.readlines() #open bed file and store it to content as list
file.close() #close file

#remove new line characters from each string in content list and split each element based on tab character
for i in range(0,len(content)):
    content[i] = content[i].replace('\n','')
    content[i] = content[i].split('\t')


occurence_list = parse_occurences(content) #scan for coverage and save info to occurence_list

#print occurence list in appropriate format to console
for i in occurence_list:
    print(str(i[0])+'\t'+str(i[1])+'\t'+str(i[2])+'\t'+str(i[3]))

exit(0)