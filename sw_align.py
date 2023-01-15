#!/usr/bin/env python3 

#Description: This program performs the Smith-Waterman local sequence alignment for 2 input fasta files and prints out the alignment with an alignment score
#             This uses a fixed scoring matrix with -1 for mismatch or gap and +1 for match

import sys

#helper function to read fasta sequence from file path
def read_sequence(seq_file_path):
    seq_file = open(seq_file_path,'r') #opens fasta file in read mode
    seq_contents = seq_file.readlines() #stores file lines in list
    seq_file.close() #closes file
    seq_string = ' ' #adds empty character at beginning (important when we align the 2D grid for alignment)
    
    #for loop scans all the lines except for the header and concatenates it to the seq_string it while removing any new line characters
    for i in range(1,len(seq_contents)):
        seq_contents[i] = str.replace(seq_contents[i],'\n','')
        seq_string += seq_contents[i]
    return seq_string #sequence string is returned

#generates a 2D grid whose dimensions is equal to the 2 sequence lengths (contains current score and previous direction (Diagonal,Left,Up))
def generate_grid(seq1,seq2):
    grid = [] #initialize grid
    for i in seq2:
        row = [] #intiialize row
        for j in seq1:
            row.append([0,'']) #append grid with defualt value of 0 and No direction
        grid.append(row)
    return grid #returns the grid

def get_cordinate_value(x,y,grid,seq1,seq2):
    x_prev = x-1
    y_prev = y-1
    coordinates = [(x_prev,y_prev),(x,y_prev),(x_prev,y)] #generates a list of the previous coordinates to scan
    vals = [(0,'')] #stores all the calculated values when adding either a match, mismatch, or gap with either of the previous grid cells, adds a defualt value of 0 and no direction since negative values not allowed
    
    #scans the list with all the previous cells (Diagonal, Left, Up)
    for i in range(0,len(coordinates)):
        if coordinates[i][0] >= 0 and coordinates[i][1] >= 0:
            #if cell exists (no negative coordinates) then proceed else ignore it
            j = i #variable J will be used in switch statement
            if seq1[x] == seq2[y] and j == 0:
                j = 3 #if you scan the diagonal (j = 0) and there is a match set J to 3
            switch = {
                0:('D',-1), #if j is 0 then its a diagnoal with mismatch
                1:('U',-1), #if j is 1 then its a gap with the previous cell being up
                2:('L',-1), #if j is 2 then its a gap with the previous cell being left
                3:('D',1) #if j is 3 then its a diagonal with a match
            }
            direction = switch.get(j) #use switch statements to determine value to add to the respective diagonal, up, or left previous cells
            vals.append((grid[coordinates[i][1]][coordinates[i][0]][0] + direction[1],direction[0])) #append the calculated score with the direction it came from to vals as a tuple
        else:
            continue
    max_score = [0,''] #set max_score with default value of 0 and No direction
    for i in range(1,len(vals)):
        #if larger max score exists update max score accordingly (thus no negative values allowed)
        if vals[i][0] >= max_score[0]:
            max_score[0] = vals[i][0]
            max_score[1] = vals[i][1]
    return max_score

#helper function to traverse through the initilaized grid and calculate the appropriate scores along with the privous cell each score was derived from       
def traverse_grid(grid,seq1,seq2):
    x = 0
    y = 0
    start = [0,0,0] #keeps track of the coordinates of the grid cell with the maximum value
    for n in grid:
        x = 0
        for m in n:
            grid[y][x] = get_cordinate_value(x, y, grid, seq1, seq2)
            if grid[y][x][0] > start[0]:
                start = [grid[y][x][0],x,y] #if grid as being traversed has a larger value update the scroe with that value and grid coordinates
            x += 1
        y+=1
    return (grid,start) #return tuple of that grid and the starting point containing the coordinates of the maximum value

#helper function to traceback the finalized grid and determine the sequence alignment along with the alignment score
def get_alignment(grid,start,seq1,seq2):

    #set the starting x and y coordinate and alignment score in the grid the same as those from the generated start list and scan from there
    y_cor = start[2]
    x_cor= start[1]
    score = start[0]
    alignment = []

    #while loop runs until you get either to the origin of the grid (0,0) coordinates or a gird cell with a value of 0 and No previous direction
    end = False
    while end == False:
        #check when to go diagonally
        if grid[y_cor][x_cor][1] == 'D':
            x_cor -= 1
            y_cor -= 1
            if seq1[x_cor+1] == seq2[y_cor+1]:
                alignment.append([seq2[y_cor+1],'|',seq1[x_cor+1]]) #if a match then add '|' chatacter
            else:
                alignment.append([seq2[y_cor+1],'*',seq1[x_cor+1]]) #if a mismatch then add a '*' character
        elif grid[y_cor][x_cor][1] == 'L':
            #check when going left
            x_cor -= 1
            alignment.append(['-',' ',seq1[x_cor+1]]) #append sequence with the gap on the other sequence
        elif grid[y_cor][x_cor][1] == 'U':
            #check when going right
            y_cor -= 1
            alignment.append([seq2[y_cor+1],' ','-']) #append sequence with the gap on the other sequence
        else:
            #break loop if origin or cell with 0 and no directon is reached (direction of previous cell direction is empty)
            end = True
    return (alignment,score) #returns the list with the alignment and the score as a tuple

#helper function to format the alignment
def format_alignment(alignment):
    seq1_format = ''
    seq2_format = ''
    align_notation = ''
    #scan from the of the alignment list towards the begining
    for i in range(len(alignment[0])-1,-1,-1):
        seq1_format += alignment[0][i][0] #stores the first sequence with the gaps
        align_notation += alignment[0][i][1] #stores the matching '|' and mismathcing '*' character accordingly
        seq2_format += alignment[0][i][2] #stores the second sequence with the gaps
    
    formated = seq1_format+'\n'+align_notation+'\n'+seq2_format #concatenates the 3 generated strings by separating them with a new line character
    return formated

#reads the sequences from 2 input fasta files from command line
seq1 = read_sequence(sys.argv[1])
seq2 = read_sequence(sys.argv[2])

#initilaize a 2D grid for the 2 sequences and travserse the grid accordingly with the correct score and driection of previous cell   
grid = generate_grid(seq1, seq2)
grid = traverse_grid(grid, seq1, seq2)
align = get_alignment(grid[0],grid[1],seq1,seq2) #traceback the grid to get the alignment
print('Score:',align[1]) #print the alignment score
print(format_alignment(align)) #print the alignment properly formated

def print_grid_space(grid,seq1,seq2):
    h_str = '     '
    for i in seq1:
       h_str += i + '          '
    print(h_str) 
    j = 0
    for i in grid:
        print(seq2[j],i)
        j+=1

#print_grid_space(grid[0], seq1, seq2)
