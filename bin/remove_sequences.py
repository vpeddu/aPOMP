import re
import sys 
import gzip
import pickle
import datetime
import subprocess
from Bio import SeqIO

lookup_dict_pickle = open(sys.argv[1], "rb")
lookup = pickle.load(lookup_dict_pickle)
print('done loading lookup table')
# needs trailing backslash
index_path = sys.argv[3]
removal_reason = str(sys.argv[4])
removal_accessions = open(sys.argv[2], 'r')
accession_list = []
# for line in removal_accessions.readlines():
#     print(line.strip())
#     accession_list.append(line.strip())
#     genus = lookup[line.strip()][1]
#     species = lookup[line.strip()][0]
#     filename = str(index_path) + str(genus) + '.genus.fasta.gz'
#     match_term = line.strip() + "|" + str(species)
#     with gzip.open(filename, "rt") as to_remove_file:
#         with open('remove_temp.fasta', 'w') as temp_fasta:
#             for record in SeqIO.parse(to_remove_file, "fasta"):
#                 if record.id != match_term:
#                     SeqIO.write(record, temp_fasta, 'fasta')
#     rename_cmd = 'mv remove_temp.fasta ' + str(genus) + '.genus.fasta'
#     zip_cmd = 'pigz ' + str(genus) + '.genus.fasta'
#     subprocess.call(rename_cmd, shell = True)
#     subprocess.call(zip_cmd, shell = True)
#     temp_fasta.close()
#     to_remove_file.close()


unique_files = {}
for line in removal_accessions.readlines():
    accession_list.append(line.strip())
    if line.strip() in lookup: 
        genus = lookup[line.strip()][1]
        if genus not in unique_files: 
            unique_files[genus] = [line.strip()]
        else: 
            unique_files[genus].append(line.strip())
for g in unique_files.keys(): 
    filename = str(index_path) + str(g) + '.genus.fasta.gz'
    with gzip.open(filename, "rt") as to_remove_file:
        with open('remove_temp.fasta', 'w') as temp_fasta:
            removal_matches = []
            for accession in unique_files[g]:
                if accession.strip() in lookup:
                    species = lookup[accession.strip()][0]
                    match_term = accession.strip() + "|" + str(species)
                    removal_matches.append(match_term)
            for record in SeqIO.parse(to_remove_file, "fasta"):
                if str(record.id) not in removal_matches:
                    SeqIO.write(record, temp_fasta, 'fasta')
                else:
                    # only wrote one
                    print('found ' + str(record.id))
                
        rename_cmd = 'mv remove_temp.fasta ' + str(genus) + '.genus.fasta'
        zip_cmd = 'pigz ' + str(genus) + '.genus.fasta'
        subprocess.call(rename_cmd, shell = True)
        subprocess.call(zip_cmd, shell = True)
        temp_fasta.close()
        to_remove_file.close()

    ts = datetime.datetime.now()
    filename = str(ts) + "_removed_accessions.txt"
    with open(filename, "w") as removal_logging_file:
        removal_logging_file.writelines('accesions removed for reason: ' + removal_reason)
        removal_logging_file.writelines('following accessions were removed: ' + str(accession_list))





# ts = datetime.datetime.now()
# filename = str(ts) + "_removed_accessions.txt"
# with open(filename, "w") as removal_logging_file:
#     removal_logging_file.writelines('accesions removed for reason: ' + removal_reason)
#     removal_logging_file.writelines('following accessions were removed: ' + str(accession_list))
