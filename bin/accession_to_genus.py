import os
import sys
import math
import time
import tqdm
import pickle
import sqlite3
import warnings
import linecache
import itertools
import multiprocessing
from taxonomy_ranks import TaxonomyRanks
from pathos.multiprocessing import ProcessingPool  

''' 
Creates a dictionary for the genus of each accession number in nucl_gb.accession2taxid from ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid

Arg 1 is nucl_gb.accession2taxid file 
Arg 2 is location of ncbi taxa.sqlite file 
'''


poolsize = 16

def _translate_merged(all_taxids, db):
    '''
    modified from ete3: https://github.com/etetoolkit/ete/blob/master/ete3/ncbi_taxonomy/ncbiquery.py
    '''
    conv_all_taxids = set((list(map(int, all_taxids))))
    cmd = 'select taxid_old, taxid_new FROM merged WHERE taxid_old IN (%s)' %','.join(map(str, all_taxids))
    #print(cmd)
    result = db.execute(cmd)
    conversion = {}
    for old, new in result.fetchall():
        conv_all_taxids.discard(int(old))
        conv_all_taxids.add(int(new))
        conversion[int(old)] = int(new)
    return conv_all_taxids, conversion

def get_lineage(taxid, db):
    '''
    modified from ete3: https://github.com/etetoolkit/ete/blob/master/ete3/ncbi_taxonomy/ncbiquery.py
    '''
    """Given a valid taxid number, return its corresponding lineage track as a
    hierarchically sorted list of parent taxids.
    
    """
    if not taxid:
        return None
    taxid = int(taxid)
    result = db.execute('SELECT track FROM species WHERE taxid=%s' %taxid)
    raw_track = result.fetchone()
    #print(raw_track)
    if not raw_track:
        #perhaps is an obsolete taxid
        _, merged_conversion = _translate_merged([taxid],db)
        if taxid in merged_conversion:
            result = db.execute('SELECT track FROM species WHERE taxid=%s' %merged_conversion[taxid])
            raw_track = result.fetchone()
        # if not raise error
        if not raw_track:
            #raw_track = ["1"]
            raise ValueError("%s taxid not found" %taxid)
        else:
            warnings.warn("taxid %s was translated into %s" %(taxid, merged_conversion[taxid]))

    track = list(map(int, raw_track[0].split(",")))
    return list(reversed(track))

def get_genus (taxid, db):
    '''
    modified from ete3: https://github.com/etetoolkit/ete/blob/master/ete3/ncbi_taxonomy/ncbiquery.py
    '''
    lin = get_lineage(taxid, db)
    query = ','.join(['"%s"' %v for v in lin])
    cmd = "select taxid, rank FROM species WHERE taxid IN (%s);" %query
    result = db.execute(cmd)
    #genuslist = []
    
    for tax, track in result.fetchall():
        if list(map(str, reversed(track.split(",")))) == ['genus']:
            #genuslist.append(tax)
            #print(genuslist)
            return tax

# use batch size to control total mem usage
chunkSize = 100

resultsDict = dict()

def chunks(lst, n):
    '''
    takes a list lst and yields list subsets of length n 
    '''
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def parseLine( line_list ):
    '''
    Takes a list of lines from nucl_gb.accession2taxid and returns the genus for each of the accessions
    '''
    # accession (with version) is second column
    ete3db = sqlite3.connect(sys.argv[2])
    #sys.stderr.write("PID:{} \n".format(os.getpid()))
    resultList = [0] * len(line_list)
    start_time = time.time()
    batchSize = len(line_list)
    for l in range(len(line_list)):
        try:
            line = line_list[l].strip()
            accession = line.split()[1]
            # taxid is third column
            tid=str(line.split()[2])
            # finding genus level taxid
            #rank_taxon = TaxonomyRanks(tid)
            #rank_taxon.get_lineage_taxids_and_taxanames()
            #genus_tid = rank_taxon.lineages[tid]['genus']
            retList = [accession, tid, get_genus(tid, ete3db)]
            # append to resultslist
            resultList[ l ] = retList
            #sys.stderr.write("PID:{} Batch completed batchsize: {} \n".format(os.getpid(), batchSize))
            #sys.stderr.write("PID:{} runtime {}  \n".format(os.getpid(),(time.time() - start_time)))
        except:
            None
            #sys.stderr.write("PID:{} FAILED".format(os.getpid()))
    if not resultList:
        print(line_list)
    return resultList

# create pool of size poolsize as defined above 
pool = ProcessingPool(nodes=poolsize)

#read nucl_gb.accession2taxid        
with open(sys.argv[1]) as nucl2gb:
    lines = nucl2gb.readlines()[1:]
print('done reading nucl2gb')

# break nucl_gb.accession2taxid into lists of size chunksize 
listoflists = list(chunks(lines, chunkSize))
print("number of batches: {}".format(len(listoflists)))
print("poolsize: {}".format(poolsize))


print('done creating chunks of size ' + str(chunkSize))

# runs genus lookup in parallel
outList = tqdm.tqdm(pool.imap( parseLine, listoflists), total = len(listoflists))

print('done finding genuses')

# create a dictionary from imap iterators where accession is key and value is [species_taxid, genus_taxid]
# wrapped in try except because not all accessions will have a valid genus
for outer in outList: 
    try:
        if outer is None: 
            print(outer)
        for inner in outer: 
            if inner[2] is not None:
                    resultsDict[inner[0]] = [inner[1],inner[2]]
    except:
        print('one failed')
print('done creating dict')
#print(resultsDict)

# write output to pickle
print('writing pickle')
with open('nucl2gb_lookup.pkl', 'wb') as handle:
    pickle.dump(resultsDict, handle, protocol=pickle.HIGHEST_PROTOCOL)