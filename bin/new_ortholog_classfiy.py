import sys 
import pysam
import taxopy
import logging
from typing import Dict, List, Tuple

taxdb = taxopy.TaxDb(nodes_dmp="nodes.dmp", names_dmp="names.dmp", merged_dmp="merged.dmp", keep_files=True)


with open(sys.argv[1]) as annotation_file:
	lines = annotation_file.readlines()

read_dict = {}
for line in lines:
    #ignore header lines
    if line.startswith('#'): 
        continue
    # get read ID 
    read_id = line.split()[0]
    # get Tax ID 
    tax_id = line.split('|')[-3].split('@')[-1]
    if read_id not in read_dict:
        read_dict[read_id] = [taxopy.Taxon(int(tax_id),taxdb)]
    else: 
        read_dict[read_id].append(taxopy.Taxon(int(tax_id),taxdb))



LOGGER = logging.getLogger()
LOGGER.setLevel(0)

class read():
    def __int__(self):
        self.id = ''
        self.cigar = ''
        self.mapq = ''
        self.seq = ''
        self.taxid = []
        self.seen = False

# acc2taxid = open(sys.argv[3], 'r')
# # TODO: append these taxids within the fasta so we don't have to build this dictionary every time
# acc_dict = {}
# for acc in acc2taxid: 
#     acc_dict[acc.split()[1]] = acc.split()[2]
# print('done building accession to taxid dict')

not_in_accs_filename = sys.argv[2] + '.accession_DNE.txt' 
not_in_accs_file = open(not_in_accs_filename, 'w')

assignments = {}
not_in_accs_file.writelines('failed LCA: \n')
for read in read_dict.keys():
    #print(len(read_dict[read].taxid))
    #print((read_dict[read].taxid))
    #print(read_dict[read].id)
    #print('\n')
    if len(read_dict[read]) <= 1:
        lca = read_dict[read][0].taxid
        #print(lca)
        if lca not in assignments:
            assignments[lca] = 1
        else: 
            assignments[lca] += 1
    else:
        # normalize alignment scores against the maximum alignment score for weighting
        #norm_ascore = [math.exp(math.sqrt(abs(i))) for i in read_dict[read].ascore]
        #norm_ascore = [(.9**(abs(i) + 1)) for i in read_dict[read].ascore]
        #norm_ascore = [float(i)/max(read_dict[read].ascore) for i in read_dict[read].ascore]
        lca_lineage = taxopy.find_majority_vote(read_dict[read], taxdb)   
        #print(lca_lineage)
        lca = lca_lineage.taxid
        #print(lca)
        if lca not in assignments:
            assignments[lca] = 1
        else: 
            assignments[lca] += 1
    # except:
    #     not_in_accs_file.write(str(read_dict[read].taxid))
            
outfilename = sys.argv[2] + '.orthologs.prekraken.tsv'

with open(outfilename, 'w') as prekraken:
    for taxa in assignments.keys():
        line = str(taxa) + '\t' + str(assignments[taxa])
        prekraken.write("%s\n" % line)
