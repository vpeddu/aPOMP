import os
import sys
import csv
import gzip
from Bio import SeqIO
from tqdm import tqdm
import multiprocessing
from ete3 import NCBITaxa
from collections import defaultdict

def parse_fasta(fasta_file):
    with gzip.open(fasta_file, "rt") as handle:
        acc_taxid = []
        genus = fasta_file.split('/')[-1].split('.')[0]
        for record in SeqIO.parse(handle, "fasta"):
            accession = record.id.split('|')[0]
            taxid = record.description.split(' ')[1]
            acc_taxid.append((record.description, taxid))
        return acc_taxid, genus

def accessions_to_taxids(accession_lookup, genus_name):
    writefile = []
    for record in accession_lookup:
        accession = record[0]
        taxid = record[1]
        #print(taxid)
        # Use taxid variable to find and save superkingdom, genus, family, species, and strain names
        try:
            lineage = ncbi.get_lineage(taxid)
        except:
            print('taxid not found:', taxid)
            continue
        ranks = ncbi.get_rank(lineage)
        names = ncbi.get_taxid_translator(lineage)
        superkingdom = None
        phylum = None
        clss = None
        order = None
        family = None
        genus = None
        species = None
        strain = None
        for rank in ranks:
            if ranks[rank] == 'superkingdom':
                superkingdom = ncbi.get_taxid_translator([rank])[rank]
            if ranks[rank] == 'phylum':
                phylum = ncbi.get_taxid_translator([rank])[rank]
            if ranks[rank] == 'class':
                clss = ncbi.get_taxid_translator([rank])[rank]
            if ranks[rank] == 'order':
                order = ncbi.get_taxid_translator([rank])[rank]
            if ranks[rank] == 'family':
                family = ncbi.get_taxid_translator([rank])[rank]
            if ranks[rank] == 'genus':
                genus = ncbi.get_taxid_translator([rank])[rank]
            if ranks[rank] == 'species':
                species = ncbi.get_taxid_translator([rank])[rank]
            if ranks[rank] == 'strain':
                strain = ncbi.get_taxid_translator([rank])[rank]
        if strain is None:
            strain = 'NaN'
    # if None not in (superkingdom, phylum, clss, order, family, genus, species):
        writefile.append([accession, superkingdom, phylum, clss, order, family, genus, species, strain])
    # else:
    #     for t in [superkingdom, phylum, clss, order, family, genus, species]:
    #         if t is None:
    #             t = 'NaN'
    #     writefile.append([accession, taxid, 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'])

    # Write the data to a CSV file
    csv_file = 'sourmash_csv_output_directory/'+ genus_name + "_sourmash_tax_table.csv"
    with open(csv_file, "w", newline="") as file:
        writer = csv.writer(file)
        writer.writerows(writefile)

def process_file(file):
    """
    Process a single file by parsing the fasta and mapping accessions to taxids.
    
    Parameters:
    - file (str): The filename to be processed.
    """
    #print('running', file)
    currentfa = fasta_dir_path + file
    t, g = parse_fasta(currentfa)
    accessions_to_taxids(t, g)

if __name__ == "__main__":
    # Check if the correct number of command line arguments are provided
    if len(sys.argv) != 2:
        print("Usage: python build-sourmash-taxtable.py <fasta_directory>")
        sys.exit(1)
    # Initialize NCBITaxa
    ncbi = NCBITaxa()
    fasta_dir_path = sys.argv[1]
    print('done reading in ncbi')
    # Read in tsv.gz file
    output_dir = "sourmash_csv_output_directory"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    fasta_dir = os.listdir(fasta_dir_path)
    with multiprocessing.Pool(processes=64) as pool:
        for _ in tqdm(pool.imap_unordered(process_file, fasta_dir, chunksize = 1), total = len(fasta_dir)):
            pass
            #for _ in tqdm.tqdm(pool.imap_unordered(do_work, tasks), total=len(tasks)):

    # for file in fasta_dir:
    #     print('running ',file)
    #     currentfa = fasta_dir_path +  file
    #     t,g = parse_fasta(currentfa)
    #     accessions_to_taxids(t,g)
    #print("Taxonomic analysis completed. Results written to", tax_table_file)
    