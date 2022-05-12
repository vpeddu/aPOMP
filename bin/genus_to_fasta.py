import os
import sys
import time
import pickle 
import subprocess
from tqdm import tqdm
from Bio import SeqIO
from ete3 import NCBITaxa


'''
Takes the accession to genus pickle (arg 1) and a folder of nt.*.tar.gz from https://ftp.ncbi.nlm.nih.gov/blast/db/ and organizes accessions into genus level fastas

Will need pigz, tar, GNU parallel, and blastcmd ainstalled locally 
'''

print('starting')
readstarttime = time.time()
lookup_dict_pickle = open(sys.argv[1], "rb")
lookup = pickle.load(lookup_dict_pickle)
print('read in lookup dictionary in ', (time.time() - readstarttime), ' seconds')
tarballs = os.listdir(sys.argv[2])

subprocess.call('mkdir genus_organized',shell=True)
# loop through each tarball in the folder of nt tars
for file in tarballs: 
    print('processing file ', file)
    #extract the tar archive 
    tar_extract_cmd = 'tar -xvzf ' + sys.argv[2] + '/' + file
    fasta_extract_cmd = 'blastdbcmd -entry all -db ' + file.split('.tar.gz')[0] + ' -out temp.fasta'
    subprocess.call(tar_extract_cmd, shell=True)
    subprocess.call(fasta_extract_cmd, shell=True)
    print('processing ' + file)
    # initialize dictionary to hold each fasta record that falls into a genus
    # key is genus, value is the seqio fasta object
    seen = {}
    with open('temp.fasta') as fastafile:
        records = SeqIO.parse(fastafile, 'fasta')
        # loop through the fasta
        for record in tqdm(records): 
            #lookup_start_time = time.time()
            # Make sure the record exists in our accession to genus lookup dictionary 
            if record.id in lookup:
                genus = lookup[record.id][1]
                record.description = lookup[record.id][0]
                #record.id = str(record.id) + "|" + str(record.description)
                # if genus is not in dictionary yet create a list containing just that record
                if genus not in seen:
                    seen[genus] = [record]
                # if genus is in dictionary append fasta record to the end of the genus fasta list
                else:
                    seen[genus].append(record) 
        print('done creating genus fasta dictionary')
        # write the genus to fasta dictionary to fastas within each genus
        for g in tqdm(seen.keys(),total = len(seen.keys())): 
            for i in seen[g]:
                i.id = str(i.id) + "|" + str(i.description)
        for g in tqdm(seen.keys(),total = len(seen.keys())): 
            tempfastafilename = 'genus_organized/' + str(g) + '.' + str(file) + '.temp.fasta'
            record_list = seen[g]
            with open(tempfastafilename, 'w') as temp_fasta:
                SeqIO.write(record_list, temp_fasta, 'fasta')
            temp_fasta.close()
        # compress all of the fastas after writing them
        compress_cmd = 'pigz -p 16 genus_organized/*.fasta'
        subprocess.call(compress_cmd, shell = True)
    # delete and close intermediates
    fasta_remove_cmd = 'rm temp.fasta'
    clean_cmd = 'rm ' + file.split('.tar.gz')[0] + '.n*'
    subprocess.call(fasta_remove_cmd, shell=True)
    subprocess.call(clean_cmd, shell=True)
    fastafile.close()
    
#concatenate all temp fastas into genus level final fastas
#print('concatenating all fastas')
#concatenate_cmd = 'ls genus_organized/ | cut -f1 -d . | sort | uniq | parallel -j 10 "cat genus_organized/{{}}.nt.* > final_index/{{}}.genus.fasta.gz"'
#subprocess.call(concatenate_cmd, shell = True)