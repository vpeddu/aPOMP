import os
import sys
import csv
import gzip
import tqdm
import time
import pickle
import subprocess
from Bio import SeqIO
from tqdm import tqdm
import multiprocessing
from ete3 import NCBITaxa
from collections import defaultdict
import IPython

# readstarttime = time.time()
# lookup_dict_pickle = open(sys.argv[1], "rb")
# lookup = pickle.load(lookup_dict_pickle)
# print('read in lookup dictionary in ', (time.time() - readstarttime), ' seconds')

def read_file_into_dict(file_path):
    lookup = {}
    with open(file_path, 'r') as file:
        for line in tqdm(file):
            columns = line.strip().split()
            key = columns[0]
            value = columns[1]
            lookup[key] = value
    return lookup

# Specify the path to the file
file_path = sys.argv[1]

# Call the function to read the file into a dictionary
lookup = read_file_into_dict(file_path)

# Create an instance of NCBITaxa
ncbi = NCBITaxa(sys.argv[2])
print('read sqlite')

tree = ncbi.get_topology([1], intermediate_nodes = True)
print('built tree, searching for genera')
genus_nodes = tree.search_nodes(rank="genus")
old_len = len(genus_nodes)

#for below: blastdbcmd -db core_nt/core_nt -entry all -outfmt %T > taxids_present.txt
taxids_file = "taxids_present.txt" 

# Read the taxids_present.txt file into a set
taxid_set = set()
with open(taxids_file, 'r') as file:
    for line in file:
        taxid = line.strip()
        taxid_set.add(int(taxid))
print('read taxids_present.txt')


def check_if_exists_already(genus_list):
    out_list = []
    for g in tqdm(genus_list, total = len(genus_list)):
        taxid = g.taxid
        fastafilename = str(taxid) + '.genus.fasta.gz'
        if not os.path.exists(os.path.join('aPOMP_genus_organized_nt', fastafilename)):
            out_list.append(g)
    return out_list


genus_nodes = check_if_exists_already(genus_nodes)
new_len = len(genus_nodes)
print(f"Removed {old_len - new_len} genera that already exist")

new_node_list = {}
for i in genus_nodes:
    new_node_list[i.taxid]=i

genera_to_keep = set()
for i in new_node_list:
    if new_node_list[i].children:
        for j in new_node_list[i].children:
            if j.rank == 'species':
                if j.taxid in taxid_set:
                    genera_to_keep.add(new_node_list[i])
print(f"{len(genera_to_keep)} genera to keep")
#IPython.embed()
genus_nodes = list(genera_to_keep.intersection(set(genus_nodes)))

# Create output directory
output_dir = "aPOMP_genus_organized_nt"
os.makedirs(output_dir, exist_ok=True)

def process_genus_node(node):
    taxid = node.taxid
    #print('running ', taxid)
    #blast_extract_cmd = f'/private/groups/kimlab/vikas/miniconda/bin/blastdbcmd -db /scratch/vpeddu/new_core_nt/core_nt -taxids {taxid} -outfmt "%f" -out {taxid}.tmp.fasta'
    fastafilename = str(taxid) + '.genus.fasta.gz'
    # if os.path.exists(os.path.join(output_dir, fastafilename)):
    #     print(f"Skipping {taxid} as it already exists")
    #     return
    try:
        fasta_file = f"{taxid}.tmp.fasta"
        #print(blast_extract_cmd)
        subprocess.run(['blastdbcmd', 
                        '-db', 
                        sys.argv[3], 
                        '-taxids', 
                        str(taxid), 
                        '-target_only',
                        '-outfmt', 
                        '"%f"', 
                        '-out', 
                        fasta_file], 
                        check=True)
        sequences = []
        with open(fasta_file, "r") as file:
            for record in SeqIO.parse(file, "fasta"):
                try:
                    species = lookup[record.id]
                    record.description = species
                    record.id = str(record.id) + "|" + str(record.description)
                    if len(record.seq) > 0:
                        sequences.append(record)
                except:
                    print(f"{record.id} not in accession lookup")
        with gzip.open(os.path.join(output_dir, fastafilename), 'wt') as temp_fasta:
            SeqIO.write(sequences, temp_fasta, 'fasta')
        os.remove(fasta_file)
    except subprocess.CalledProcessError as e:
        #print(f" {taxid} failed with return code {e.returncode}")
        #os.remove(fasta_file) UNCOMMENT WHEN DONE TESTING
        with open("failed_taxa.txt", "a") as file:
            file.write(str(taxid) + "\n")       
        os.remove(fasta_file)
# Create a multiprocessing pool
pool = multiprocessing.Pool(24)

# Map the process_genus_node function to each genus node in parallel
results = []
for result in tqdm(pool.imap_unordered(process_genus_node, genus_nodes, chunksize=5), total=len(genus_nodes)):
    results.append(result)

# Close the pool to prevent any more tasks from being submitted
pool.close()

# Wait for all processes to complete
pool.join()


