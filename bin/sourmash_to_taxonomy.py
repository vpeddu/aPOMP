# script to take genus names from sourmash lca summarize and return a file with list of genus level taxids

import csv
import sys 
import pandas as pd

print(sys.version)
from ete3 import NCBITaxa

ncbi = NCBITaxa(dbfile = 'taxa.sqlite')

with open(sys.argv[1], 'r') as file:
    reader = csv.reader(file)
    genus_taxids = []
    for row in reader:
        # column 6 is genus
        # first row is header so skipping it here
        if row[6] and row[6] != 'genus':
            #print(row[6], ncbi.get_name_translator[str(row[6])])
            # not all taxa are up to date in the NCBI sqlite database. Small number but needs to be fixed 
            tid = ncbi.get_name_translator([str(row[6])])
            if tid: 
                genus_taxids.append(tid[row[6]][0])

df = pd.DataFrame(genus_taxids) 
filename = sys.argv[2] + '.sourmash_to_genus.txt'
df.to_csv(filename, index = False) 

