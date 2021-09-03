import sys
import pickle
from Bio import SeqIO

with open(sys.argv[1]) as nucl2gb:
	accs = nucl2gb.readlines()
print('done reading nucl2gb')

lookup = {}

for acc in accs:
	if acc.split()[2] not in lookup:
		lookup[acc.split()[2]] = [acc.split()[1]]
	else: 
		lookup[acc.split()[2]].append(acc.split()[1])

with open('nucl2gb_lookup.pkl', 'wb') as handle:
    pickle.dump(lookup, handle, protocol=pickle.HIGHEST_PROTOCOL)

print('done creating accession to taxid lookup')

record_dict = SeqIO.to_dict(SeqIO.parse(sys.argv[2], "fasta"))
#print(record_dict['NC_056724.1'].description)  # use any record ID
with open('ref_nt.pkl', 'wb') as handle:
    pickle.dump(record_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
print('done creating fasta pkl')
