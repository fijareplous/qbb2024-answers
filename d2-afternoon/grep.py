#!/usr/bin/env python3

import sys

file = open(sys.argv[1])

for line in file:
    line = line.rstrip("\n")
    print(line)

file.close()




