import sys
import math
import pysam
import numpy
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
        self.alen = ''
        self.qlen = ''
        self.weights = ''


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
        #read_dict[record.query_name].seq = record.query_sequence
        # store taxid in a list as taxopy object
        # this is necessary for LCA later on 
        read_dict[record.query_name].taxid = [record_tid]
        # read has been seen (probably don't need this)
        read_dict[record.query_name].seen = True
        # number of bases overlapping with reference
        read_dict[record.query_name].alen = [record.query_alignment_length]
        read_dict[record.query_name].qlen = record.query_length
    # if read aready exists in read dictionary 
    else: 
        #if not read_dict[record.query_name].seq: 
            # read_dict[record.query_name].seq = record.query_sequence
        # append alignment scores
        read_dict[record.query_name].ascore.append(record.get_tag("AS"))
        # store mapq
        read_dict[record.query_name].mapq.append(record.query_qualities)
        # append this taxid to the read taxid list as taxopy object
        read_dict[record.query_name].taxid.append(record_tid)
        # 
        read_dict[record.query_name].alen.append(record.query_alignment_length)
print('done creating read dictionary')

for read in read_dict.keys():
    if len(read_dict[read].alen) > 1:
        indexed_overlap_sort = numpy.argsort(read_dict[read].alen)
        read_dict[read].alen = numpy.array(read_dict[read].alen)
        top_10 = indexed_overlap_sort[::-1][0:10]
        taxid_list = numpy.array(read_dict[read].taxid)
        #print(taxid_list[indexed_overlap_sort][::-1][0:10])
        read_dict[read].taxid = taxid_list[top_10]
        read_dict[read].weights = top_10
        #print(read_dict[read].alen,read_dict[read].alen[numpy.array([indexed_overlap_sort])][::-1])


def most_frequent(List):
    return max(set(List), key = List.count)
    
assignments = {}
not_in_accs_file.writelines('failed LCA: \n')
for read in read_dict.keys():
    #print(len(read_dict[read].taxid))
    #print((read_dict[read].taxid))
    #print(read_dict[read].id)
    #print('\n')
    taxopy_read_list = []
    for tid in read_dict[read].taxid: 
        taxopy_read_list.append(taxopy.Taxon(int(tid),taxdb))
    if len(taxopy_read_list) <= 1:
        lca = taxopy_read_list[0].taxid
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
        
        #if taxopy_read_list.count(most_frequent(taxopy_read_list)) / len(taxopy_read_list) > .5:
        #        lca = most_frequent(taxopy_read_list)
        #else:
        #        lca_lineage = taxopy.find_majority_vote(taxopy_read_list, taxdb)   
        lca_lineage = taxopy.find_majority_vote(taxopy_read_list, taxdb, weights = read_dict[read].weights.tolist())   
        #print(lca_lineage)
        lca = lca_lineage.taxid
        # if lca == 543:
        #     print(taxopy_read_list)
        #     print('\n')
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
#from IPython import embed; embed()