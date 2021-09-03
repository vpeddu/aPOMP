import sys
from Bio import SeqIO
from ete3 import NCBITaxa

ncbi = NCBITaxa()

with open("nucl2gb_lookup.pkl") as nucl2gb:
	accs = nucl2gb.readlines()
print('done reading nucl2gb')

lookup = {}
for acc in accs:
	if acc.split()[2] not in lookup:
		lookup[acc.split()[2]] = [acc.split()[1]]
	else: 
		lookup[acc.split()[2]].append(acc.split()[1])
print('done creating accession to taxid lookup')

record_dict = SeqIO.to_dict(SeqIO.parse("ref_nt.pkl", "fasta"))
#print(record_dict['NC_056724.1'].description)  # use any record ID

print('done creating fasta lookup')

with open(sys.argv[1]) as krakenreport:
	lines = krakenreport.readlines()

genus_list = []
for i in lines: 
	if i.split()[3] == 'G':
		genus_list.append(i.split()[4])

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