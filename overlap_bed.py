#!/usr/bin/env python3

#Description: This program uses 2 BED files as input and then tries to find which of the regions overlap by at least a certain percentage threshold set by the user. 
#             Optimiziations such as the BED files being sorted first are used to speed up the lookup time.

import argparse

#processes arguments such as the input files, output file, match percentage, and join
parser = argparse.ArgumentParser(description="Prints overlapping regions from 2 bed files that have a minimum overlap percentage set by the user")
parser.add_argument('-i1',metavar='--input1',type=str,help='Enter the first BED file to detect overlaps with the second one')
parser.add_argument('-i2',metavar='--input2',type=str,help='Enter the second BED file')
parser.add_argument('-m',metavar='--overlap',type=int,help='Enter the minimum overlap percentage 0-100')
parser.add_argument('-o',metavar='--output',type=str,help='Enter the output BED file name containing the overlaps of the first BED file contrasted to the second one')
parser.add_argument('-j',action='store_true',help='Join matching lines from both bed files') #flag -j does not take any values

args = parser.parse_args()

#helper function to check if arguments are valid by passing a list. If one of the values is false program terminates
def check_args(arg_list):
    for i in arg_list:
        if i == False:
            exit(1)
args_list = [True,True,True] #list that checks if input and match arguments are valud

#helper function to open the bed files and return a dictionary with chromosome names as key and list of all the start-end postions as value
def open_bed(path):
    file = open(path,'r') #opens filepath given in read mode
    lines = file.readlines() #stores each line of the bed file into a list
    chr_list = {} #initialize dictionary to store chromosom and start-end postion info
    file.close() #close bed file

    #for loop pairs the corresponding chromosome names with its associated start-end position pairs and adds that to the dictionary
    for i in range(0,len(lines)):
        lines[i] = lines[i].replace('\n','')#removed the new line character
        lines[i] = lines[i].split('\t')#slpits each line based on the tab delimeter

        #if chromosome name not already a key in dicitonary add it with an empty list as value
        if lines[i][0] not in chr_list:
            chr_list[lines[i][0]] = []
        
        #convert each start and end postions to integers
        lines[i][1] = int(lines[i][1])
        lines[i][2] = int(lines[i][2])
        #add each start-end pair along with additional info for that line as a list and append that list to the dictionary
        chr_list[lines[i][0]].append(lines[i][1:])

    #for loop appropriately sorts the start-end postions among the dictionary (important for optimiziation)
    for i in chr_list.keys():
        chr_list[i].sort()
    return chr_list #returns the dictionary with the needed info for each chromosome

#helper function that detects overlaps of the first input dictionary from the first bed file when contrasted to the second one that passes a specified matching percentage threshold
#returns list of overlaped positions
def detect_overlaps(input1, input2, match):
    out_list = [] #list of overlaped positions
    

    #traverse through each key in dictionary derived from first bed file
    for i in input1.keys():
        #Each start-end positions in dictionary are sorted (plays big role in optimizing search time). A begin_index variable is used to keep track of the postion of the last match.
        #since new start-end postion from input1 will have a larger value it will only need to start scanning at the matching postion of input2 for the previous start-end postion
        #instead at the begining of the input2. Also if start postion of input2 is larger than end position of input1 the search stops since the next start-end postions of input2
        #will be even larger hence never have an overlap with input1.

        begin_index = 0 #begin_index is initialized to 0 when working with new chromosome

        #traverse through each start-end pair for each chromosome
        for j in input1[i]:
            found = False #becomes true if overlap passing matching percentage threshold is found, intially false

            #find overlaps for current start-end position of input1 in input2, start at begin_index since previous start-end postions in input2 are smaller hence no overlap
            for k in range(begin_index,len(input2[i])):

                #if current end position from input1 is larger than start position from input2 then exit loop as no new matches will be found
                if j[1] < input2[i][k][0]:
                    #if no overlap has been found for start-end position from input1 update begin_index
                    if found == False:
                        begin_index = k
                    break #exit loop

                overlap = min(j[1],input2[i][k][1]) - max(j[0],input2[i][k][0]) #caluclate bases overlapping (should be larger than 0 if there is some overlap)
                length = j[1] - j[0] #calculate length of chromosome segment from input1

                #if length is 0 exit the loop as no overlap will ever be found
                if length == 0:
                    break

                percentage = (float(overlap)/float(length))*100 #calculate the percentage of segment from input1 overlaping in input2

                #if some overlap exists
                if overlap > 0:
                    #if it is first overlap then set found to True and update begin_index
                    if found == False:
                        begin_index = k
                        found = True
                    #if the percentage passes the matching threshold
                    if percentage >= match:
                        #if join flag is present append the the list of overlapped postions the chromosme name, start-end position from input1 and where it overlaps in input2
                        if args.j:
                            out_list.append([i,j,input2[i][k]])
                        else:
                            #if join flag is absent check if the current start-end postion for that chromosome already is in list
                            if len(out_list)>0 and out_list[len(out_list)-1][0] == i and out_list[len(out_list)-1][1][0] == j[0] and out_list[len(out_list)-1][1][1] == j[1]:
                                pass #if in list do not add it again as it is redundant
                            else:
                                #if not in list add the chromosme name and the start-end position from input1
                                out_list.append([i,j])
    return out_list #return the overlaped list

#helper function to generate the output containing overlaped info in correct format and given a file name and overlap list and a README file with the number of overlapped matches
def generate_output_file(out_name,out_content):
    out_file = open(out_name,'w') #create the output file with write permssions
    if args.j:
        #if join flag is present
        for i in out_content:
            #for each sub-list in the overlap list
            l1 = '' #will include info from overlapped postions from first bed file
            l2 = '' #will include info from second bed file that first bed file had overlap
            #tab delimit line with start-end postion info related to first bed file 
            for j in i[1]:
                l1 += str(j) + '\t'
            #tab delimit line with start-end postion info related to second bed file 
            for j in i[2]:
                l2 += str(j) + '\t'     
            #write line to output file with chromosome name, start-end from first bed file, and start-end from second bed file with overlap
            out_file.write(str(i[0])+'\t'+l1+str(i[0])+'\t'+l2+'\n')    
    else:
        #if join flag is absent
        for i in out_content:
            #for each sub-list in the overlap list
            l1 = '' #will include info from overlapped postions from first bed file
            #tab delimit line with start-end postion info related to first bed file 
            for j in i[1]:
                l1 += str(j) + '\t'  
            #write line to output file with chromosome name, start-end from first bed file woth overlap
            out_file.write(str(i[0])+'\t'+l1+'\n')       
    out_file.close() #close output file once done
    readme = open('README.txt','w') #create a readme file with write access
    #write number of overlaped matches from first bed file (if join flag present, duplicates with overlap in mutliple regions to second bed file will also be included)
    readme.write('Number of segments from '+str(args.i1)+' file overlaping in '+str(args.i2)+' file\nwith a marching percent threshold of '+str(args.m)+' percent is '+str(len(out_content)))
    readme.close() #close readme file



#try to open the first bed file, process its contents and store them in a dictionary, if it cannot open it an expection is thrown
try:
   input1 = {} 
   input1 = open_bed(args.i1)
except IOError as x:
    print('Could not open the first BED file named '+str(args.i1)+'!')
    args_list[0] = False

#try to open the second bed file, process its contents and store them in a dictionary, if it cannot open it an expection is thrown
try:
    input2 = {}
    input2 = open_bed(args.i2)
except IOError as x:
    print('Could not open the second BED file named '+str(args.i2)+'!')
    args_list[1] = False

#checks if matching percentage is in a valid range
if args.m < 0 or args.m > 100:
    print('Match percentage is invalid. Must be between 0 and 100!')
    args_list[2] = False

#checks if input file and matching percentage arguments are valid if not program terminates
check_args(args_list)

out = detect_overlaps(input1,input2,args.m) #detect overlaps of first bed file to the second one and return result to out

generate_output_file(args.o,out) #generate ouput file and readme file in proper format based on contents of out with the provided output file name

exit(0) #program exists succesfully with no errors
