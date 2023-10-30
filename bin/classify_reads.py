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
        self.refname = ''
        self.reflen = ''
        

def ref_longer_than_read(al, ql ): 
    proportion_read_used =  al / ql
    return proportion_read_used

def read_longer_than_ref(al, rl):  
    proportion_read_used = al / rl
    return proportion_read_used

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
        #read_dict[record.query_name].mapq = [record.query_qualities]
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
        read_dict[record.query_name].qlen = record.infer_query_length()
        read_dict[record.query_name].refname = [record.reference_name]
        read_dict[record.query_name].reflen = [record.reference_length]
    # if read aready exists in read dictionary
    else:
        # append alignment scores
        read_dict[record.query_name].ascore.append(record.get_tag("AS"))
        # store mapq
        #read_dict[record.query_name].mapq.append(record.query_qualities)
        # append this taxid to the read taxid list as taxopy object
        # this taxid is plasmid
        if acc_num.startswith('p_'):
            read_dict[record.query_name].taxid.append(36549)
        else:
            read_dict[record.query_name].taxid.append(record_tid)
        # append length of alignment
        read_dict[record.query_name].alen.append(record.query_alignment_length)
        
        read_dict[record.query_name].refname.append(record.reference_name)
        
        read_dict[record.query_name].reflen.append(record.reference_length)
print('done creating read dictionary')

# for each read, if there is more than one hit per read, weight the top 10 alignments by the length of their aligned sequence
# used later on to weight the LCA to give spurious alignments lower priority
for read in read_dict.keys():
    if True:
        weights_list = []
        if len(read_dict[read].alen) > 1: # and '36549' not in read_dict[read].taxid:
            for aln in range(len(read_dict[read].reflen)): # for each hit to a read
                #print(read_dict[read].qlen)
                #print(read_dict[read].reflen[aln])
                #print(aln)
                if int(read_dict[read].reflen[aln] / read_dict[read].qlen  > 0.1): # if the length of the reference is not 40% of the length of the read, assign the lowest weight
                    if int(read_dict[read].qlen) > int(read_dict[read].reflen[aln]):
                        #weights_list.append(read_longer_than_ref(read_dict[read].alen[aln], read_dict[read].qlen))
                        weights_list.append(read_dict[read].ascore[aln])
                    else:
                        #weights_list.append(ref_longer_than_read(read_dict[read].alen[aln], read_dict[read].reflen[aln]))
                        weights_list.append(read_dict[read].ascore[aln])
                else: # if the length of the reference is not 40% of the length of the read, assign the lowest weight
                    weights_list.append(0)
            #print(weights_list)
            #print(read_dict[read].taxid)
            #print(read_dict[read].id)
            indexed_overlap_sort = numpy.argsort(weights_list) # get sort positions
            indexed_overlap_sort = numpy.array(indexed_overlap_sort) # create numpy array of aligned lengths
            top_10 = indexed_overlap_sort[::-1][0:10] # order the top 10 aligned lengths backwards
            taxid_list = numpy.array(read_dict[read].taxid) #get all taxids
            read_dict[read].taxid = taxid_list[top_10] # extract the taxids for the top sorted aligned lengths
            read_dict[read].weights = numpy.array([read_dict[read].ascore[pos] for pos in top_10])
            #old read_dict[read].weights = top_10 # assign taxids to the weights


    else: 
        if len(read_dict[read].alen) > 1:
            indexed_overlap_sort = numpy.argsort(read_dict[read].alen) # get sort positions
            read_dict[read].alen = numpy.array(read_dict[read].alen) # create numpy array of aligned lengths
            top_10 = indexed_overlap_sort[::-1][0:10] # order the top 10 aligned lengths backwards 
            taxid_list = numpy.array(read_dict[read].taxid) #get all taxids
            read_dict[read].taxid = taxid_list[top_10] # extract the taxids for the top sorted aligned lengths
            read_dict[read].weights = top_10 # assign taxids to the weights

# reweights hits to strains more heavily than spcies to push towards strain level specificity
def weight_strains(r,treadlist, weight_multiplier,ll):
    straincount = {}
    for t in range(len(treadlist)):
        if treadlist[t].rank == 'strain':
            #tax = taxopy.taxid_from_name(treadlist[t].name,taxdb)[0]
            #if ll.name in [x for x in treadlist[t].name_lineage]:
                # if tax not in straincount:
                #     straincount[tax] = 1
                # else:
                #     straincount[tax] += 1
                #     print(r,read_dict[r].refname, ll.name)
            read_dict[r].weights[t] = read_dict[r].weights[t] * weight_multiplier

# look into collections.defaultdict
assignments = {}
taxid_to_read = {}
read_id_to_taxid = {}
#not_in_accs_file.writelines('failed LCA: \n')
for read in read_dict.keys():
    #print(read)
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
        if numpy.any(read_dict[read].weights <= 0): # if any of the alignment scores are negative or 0 
            if numpy.all(read_dict[read].weights <= 0):  # check if all alignment scores are negative or 0
                lca = 0 # assign the read as unclassified 
            else:
                read_dict[read].weights[read_dict[read].weights < 0] = 0 # replace only the negative alignment scores with weight 0
        else:
            lca_lineage = taxopy.find_majority_vote(taxopy_read_list, taxdb, weights = read_dict[read].weights.tolist())
            #weight_strains(read, taxopy_read_list, 3, lca_lineage)
            lca = lca_lineage.taxid
            #old lca_lineage = taxopy.find_majority_vote(taxopy_read_list, taxdb, weights = read_dict[read].weights.tolist())
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
    for taxa in taxid_to_read.keys():
        dirname = str(taxa) + '.reads'
        if not os.path.exists(dirname):
            os.makedirs(dirname)
        tmp_singlebam_filename = dirname + '/' + str(taxa) + '.read_ids.txt'
        with open(tmp_singlebam_filename, 'w') as read_id_write:
            for line in taxid_to_read[taxa]:
                read_id_write.write(f"{line}\n")

#    print(assignments)
#    print(taxid_to_read)

# For testing
#from IPython import embed; embed()