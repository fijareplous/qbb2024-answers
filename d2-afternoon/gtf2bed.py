#!/usr/bin/env python3

import sys

test = ["a", "b", "c", "d", "e"]

my_file = open(sys.argv[1])

for line in my_file:
    line = line.strip().split()
    print((line[0] + "\t" + line[3] + "\t" + line[4] + "\t" + line[13].replace('"','').replace(";",'')))

my_file.close()

