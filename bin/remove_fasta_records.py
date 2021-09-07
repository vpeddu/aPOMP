#!/usr/bin/env python
# -*- coding: utf-8 -*-
# from https://www.biostars.org/p/4881/


import sys
from Bio import SeqIO

fasta_file = sys.argv[1]  # Input fasta file
remove_file = sys.argv[2] # Input wanted file, one gene name per line
result_file = sys.argv[3] # Output fasta file

remove = set()
with open(remove_file) as f:
    for line in f:
        line = line.strip()
        if line != "":
            remove.add(line)

fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')

with open(result_file, "w") as f:
    for seq in fasta_sequences:
        nuc = seq.seq.tostring()
        if nuc not in remove and len(nuc) > 0:
            SeqIO.write([seq], f, "fasta")