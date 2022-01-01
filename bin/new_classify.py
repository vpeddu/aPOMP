import sys
import math
import pysam
import taxopy
import logging

LOGGER = logging.getLogger()
LOGGER.setLevel(0)

class read():
    def __int__(self):
        self.id = ''
        self.ascore = ''
        self.mapq = ''
        self.seq = ''
        self.taxid = []
        self.seen = False


#print('done building accession to taxid dict')

bamfile = pysam.AlignmentFile(sys.argv[1], "rb")
print('done reading in bamfile')
read_dict = {}

not_in_accs_filename = sys.argv[2] + '.accession_DNE.txt' 
not_in_accs_file = open(not_in_accs_filename, 'w')

# load taxdb
#taxdb = taxopy.TaxDb(nodes_dmp="nodes.dmp", names_dmp="names.dmp")
taxdb = taxopy.TaxDb(nodes_dmp="nodes.dmp", names_dmp="names.dmp", merged_dmp="merged.dmp", keep_files=True)
# Loop through reach record in the bam and store info from each read
# Same read may appear multiple times in a bam since it can hit multiple references
for record in bamfile: 
    # accession number is first split in the reference name
    acc_num = record.reference_name.split('|')[0]
    #print(record.reference_name)
    # Taxid is second split in reference name
    record_tid = record.reference_name.split('|')[1]
    #print(record_tid)
    # if read ID has not been seen before
    if record.query_name not in read_dict:
        # create read object
        read_dict[record.query_name] = read()
        # store read in read dictionary
        read_dict[record.query_name].id = record.query_name
        # store alignment score
        # used to weight LCA below
        read_dict[record.query_name].ascore = [record.get_tag("AS")]
        # store mapq
        read_dict[record.query_name].mapq = [record.query_qualities]
        # store sequence (probably don't need this)
        read_dict[record.query_name].seq = record.query_sequence
        # store taxid in a list as taxopy object
        # this is necessary for LCA later on 
        read_dict[record.query_name].taxid = [taxopy.Taxon(int(record_tid),taxdb)]
        # read has been seen (probably don't need this)
        read_dict[record.query_name].seen = True
    # if read aready exists in read dictionary 
    else: 
        # append alignment scores
        read_dict[record.query_name].ascore.append(record.get_tag("AS"))
        # store mapq
        read_dict[record.query_name].mapq.append(record.query_qualities)
        # append this taxid to the read taxid list as taxopy object
        read_dict[record.query_name].taxid.append(taxopy.Taxon(int(record_tid),taxdb))


assignments = {}
not_in_accs_file.writelines('failed LCA: \n')
for read in read_dict.keys():
    #print(len(read_dict[read].taxid))
    #print((read_dict[read].taxid))
    #print(read_dict[read].id)
    #print('\n')
    if len(read_dict[read].taxid) <= 1:
        lca = read_dict[read].taxid[0].taxid
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
        lca_lineage = taxopy.find_majority_vote(read_dict[read].taxid, taxdb, weights=norm_ascore,)   
        #print(lca_lineage)
        lca = lca_lineage.taxid
        #print(lca)
        if lca not in assignments:
            assignments[lca] = 1
        else: 
            assignments[lca] += 1
    # except:
    #     not_in_accs_file.write(str(read_dict[read].taxid))
            
outfilename = sys.argv[2] + '.prekraken.tsv'

with open(outfilename, 'w') as prekraken:
    for taxa in assignments.keys():
        line = str(taxa) + '\t' + str(assignments[taxa])
        prekraken.write("%s\n" % line)
