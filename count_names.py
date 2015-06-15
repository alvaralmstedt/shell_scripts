#!/usr/bin/env python

import sys

infile = open(sys.argv[1], "r")
counts = {}

with infile:
	for line in infile.readlines():
		line = line.rstrip("\n")
		if line in counts:
			counts[line] += 1
		else:
			counts[line] = 1

for name in counts:
	print name, counts[name]
