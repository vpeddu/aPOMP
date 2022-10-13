import sys 
import csv

with open(sys.argv[1]) as f:
    lines = f.readlines()
ID_counts = {}
for i in lines:
    clean_up = i.strip()
    separate = clean_up.split()
    if separate[0] in ID_counts: 
        ID_counts[separate[0]] = int(separate[1]) + int(ID_counts[separate[0]])
    else:
        ID_counts[separate[0]] = separate[1]

outfilename = str(sys.argv[2]) + '.merged.prekraken.tsv'

with open(outfilename, 'w') as csv_file:  
    writer = csv.writer(csv_file, delimiter = '\t')
    for key, value in ID_counts.items():
        writer.writerow([key, value])

#with open(outfilename, 'w') as x:
#    x.write(product) 