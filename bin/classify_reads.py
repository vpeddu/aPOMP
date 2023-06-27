import os
import sys
import math
import pysam
import numpy
import string
import random
import taxopy
import logging

LOGGER = logging.getLogger()
LOGGER.setLevel(0)

# class to hold data for each read
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


# read in bamfile
bamfile = pysam.AlignmentFile(sys.argv[1], "rb")
print('done reading in bamfile')

# initialize dictionary to hold each read object
read_dict = {}

not_in_accs_filename = sys.argv[2] + '.accession_DNE.txt'
#not_in_accs_file = open(not_in_accs_filename, 'w')

# load taxdb
#taxdb = taxopy.TaxDb(nodes_dmp="nodes.dmp", names_dmp="names.dmp")
taxdb = taxopy.TaxDb(nodes_dmp="nodes.dmp", names_dmp="names.dmp", merged_dmp="merged.dmp", keep_files=True)
# Loop through reach record in the bam and store info from each read
# Same read may appear multiple times in a bam since it can hit multiple references
for record in bamfile:
    # accession number is first split in the reference name
    acc_num = record.reference_name.split('|')[0]
    # Taxid is second split in reference name
    record_tid = record.reference_name.split('|')[1]
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
        # store taxid in a list as taxopy object
        # this is necessary for LCA later on
        # if the record hits a plasmid, assign it the plasmid taxid 36549 not the species taxid
        if acc_num.startswith('p_'):
            read_dict[record.query_name].taxid = [36549]
        else:
            read_dict[record.query_name].taxid = [record_tid]
        # read has been seen (probably don't need this)
        read_dict[record.query_name].seen = True
        # number of bases overlapping with reference
        read_dict[record.query_name].alen = [record.query_alignment_length]
        read_dict[record.query_name].qlen = record.query_length
    # if read aready exists in read dictionary
    else:
        # append alignment scores
        read_dict[record.query_name].ascore.append(record.get_tag("AS"))
        # store mapq
        read_dict[record.query_name].mapq.append(record.query_qualities)
        # append this taxid to the read taxid list as taxopy object
        # this taxid is plasmid
        if acc_num.startswith('p_'):
            read_dict[record.query_name].taxid.append(36549)
        else:
            read_dict[record.query_name].taxid.append(record_tid)
        # append length of alignment
        read_dict[record.query_name].alen.append(record.query_alignment_length)
print('done creating read dictionary')

# for each read, if there is more than one hit per read, weight the top 10 alignments by the length of their aligned sequence
# used later on to weight the LCA to give spurious alignments lower priority
for read in read_dict.keys():
    if len(read_dict[read].alen) > 1:
        indexed_overlap_sort = numpy.argsort(read_dict[read].alen)
        read_dict[read].alen = numpy.array(read_dict[read].alen)
        top_10 = indexed_overlap_sort[::-1][0:10]
        taxid_list = numpy.array(read_dict[read].taxid)
        read_dict[read].taxid = taxid_list[top_10]
        read_dict[read].weights = top_10


def most_frequent(List):
    return max(set(List), key = List.count)

# look into collections.defaultdict
assignments = {}
taxid_to_read = {}
read_id_to_taxid = {}
#not_in_accs_file.writelines('failed LCA: \n')
for read in read_dict.keys():
    taxopy_read_list = []
    for tid in read_dict[read].taxid:
        taxopy_read_list.append(taxopy.Taxon(int(tid),taxdb))
    # if there is only one hit for a read
    if len(taxopy_read_list) <= 1:
        lca = taxopy_read_list[0].taxid
        #print(lca)
        if lca not in assignments:
            assignments[lca] = 1
        else:
            assignments[lca] += 1
        if lca not in taxid_to_read:
            taxid_to_read[lca] = [read_dict[read].id]
        else:
            taxid_to_read[lca].append(read_dict[read].id)
        if sys.argv[3] == 'save':
            read_id_to_taxid[read_dict[read].id] = lca
    # if there are multiple hits for a read, run taxopy LCA weighted by alignment lengths
    else:
        lca_lineage = taxopy.find_majority_vote(taxopy_read_list, taxdb, weights = read_dict[read].weights.tolist())
        lca = lca_lineage.taxid
        if lca not in assignments:
            assignments[lca] = 1
        else:
            assignments[lca] += 1
        if lca not in taxid_to_read:
            taxid_to_read[lca] = [read_dict[read].id]
        else:
            taxid_to_read[lca].append(read_dict[read].id)
        if sys.argv[3] == 'save':
            read_id_to_taxid[read_dict[read].id] = lca
# write the number of hits per taxa to this temp file for krakenuniq-report to turn into a pavian report later on
outfilename = sys.argv[2] + '.prekraken.tsv'

with open(outfilename, 'w') as prekraken:
    for taxa in assignments.keys():
        line = str(taxa) + '\t' + str(assignments[taxa])
        prekraken.write("%s\n" % line)

with open('taxid_to_read.csv', 'w') as prekraken:
    for taxa in taxid_to_read.keys():
        line = str(taxa) + '\t' + str(taxid_to_read[taxa])
        prekraken.write("%s\n" % line)

if sys.argv[3] == 'save':
	bamfile = pysam.AlignmentFile(sys.argv[1], "rb")
	for record in bamfile:
		if record.query_name in read_id_to_taxid:
			if not os.path.exists(str(read_id_to_taxid[record.query_name])):
				os.makedirs(str(read_id_to_taxid[record.query_name]))
			lca = str(read_id_to_taxid[record.query_name])
			tmp_singlebam_filename = lca + '/' + lca + '.' +  str(record.query_name) + '.singlebam.bam'
			tmp_sb_out = pysam.AlignmentFile(tmp_singlebam_filename.join(random.choices(string.ascii_letters, k=10)), 
                                    template = bamfile,
                                    mode= 'wb')
			tmp_sb_out.write(record)
			tmp_sb_out.close()
bamfile.close()
#    print(assignments)
#    print(taxid_to_read)

# For testing
#from IPython import embed; embed()