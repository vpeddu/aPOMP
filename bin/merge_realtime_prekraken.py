import sys 

with open(sys.argv[1]) as f:
    lines = f.readlines()
print(lines)
ID_counts = {}
for i in lines:
    clean_up = i.strip()
    separate = clean_up.split()
    if separate[0] in ID_counts: 
        ID_counts[separate[0]] = int(separate[1]) + int(ID_counts[separate[0]])
    else:
        ID_counts[separate[0]] = separate[1]
product = str(ID_counts)

outfilename = str(sys.argv[2]) + '.merged.prekraken.tsv'
with open(outfilename, 'w') as x:
    x.write(product)