#!/usr/bin/env python3
#Grep is like a search function; prints lines in a file where it matches

import sys

file = open(sys.argv[2])
#Allows us to type in any file into terminal, sys.argv is an attribute/list, the list contains 
#instructions, [./grep.py, '<word being searched>', '<file name>'] 

pattern = sys.argv[1]
#Specifies which word/phrase we're trying to search for

for line in file:
    line = line.rstrip("\n")
    #If pattern exists in the line:
    if pattern in line: 
        print(line)

file.close()



