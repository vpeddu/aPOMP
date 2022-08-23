import sys
import pickle5 as pickle
from Bio import SeqIO
from ete3 import NCBITaxa

ncbi = NCBITaxa()
nucl_pkl = sys.argv[2] + "/nucl2gb_lookup.pkl"
with open(nucl_pkl, 'rb') as nucl2gb:
	lookup = pickle.load(nucl2gb)
print('done reading nucl2gb')

# lookup = {}
# for acc in accs:
# 	if acc.split()[2] not in lookup:
# 		lookup[acc.split()[2]] = [acc.split()[1]]
# 	else: 
# 		lookup[acc.split()[2]].append(acc.split()[1])
# print('done creating accession to taxid lookup')

refnt_pkl = sys.argv[2] + "/ref_nt.pkl"
#record_dict = SeqIO.to_dict(SeqIO.parse(refnt_pkl, "fasta"))
with open(refnt_pkl, 'rb') as refnt:
    record_dict = pickle.load(refnt)
#print(record_dict['NC_056724.1'].description)  # use any record ID

print('done creating fasta lookup')

with open(sys.argv[1]) as krakenreport:
	lines = krakenreport.readlines()

genus_list = []

genus_count  = 0

for i in lines: 
	if i.split()[3] == 'G':
		genus_list.append(i.split()[4])
		genus_count = genus_count + 1

if genus_count > 0: 
    print('found ' + str(genus_count) +' genera')
else: 
    print('no genera found. exiting')
    sys.exit(1)

print('genus list:')
print(genus_list)

descendants = []
for genus in genus_list: 
	 descendants += (ncbi.get_descendant_taxa(genus,intermediate_nodes=True))
#print(descendants)
#print(descendants)
# print(31754 in descendants)
with open('species.fasta', 'w') as newfasta:
	for taxa in descendants: 
		try:
			records = lookup[str(taxa)]
			for fasta in records:
				record_dict[fasta].id += '|' + str(taxa)
				#records_dict[fasta].description = ''
				SeqIO.write(record_dict[fasta], newfasta, 'fasta')
		except:
			None	

newfasta.close()