import os
import re
import sys
import csv
import taxoniq
import numpy as np
import pandas as pd
from IPython.display import display



class sample(): 
    def __init__(self, name, species_list):
        self.name = name
        self.unclassified = 0
        self.taxon_list = self.populate_lineage(species_list)
        self.correctly_classified = {'order':[],'family':[], 'genus':[], 'species':[]}
        self.total_classified_count = 0
        self.precision = 0
        self.recall = 0
        self.f1 = 0

    def populate_lineage(self, species_list):
        taxon_dict = {'order':[],'family':[], 'genus':[], 'species':[]}
        for spc in species_list: 
            current_species = taxoniq.Taxon(spc)
            lineage = current_species.ranked_lineage
            for taxa in lineage:
                rank =  str(taxa.rank.name)
                if rank in taxon_dict: 
                    taxon = int(str(taxa).split('(')[1].split(')')[0])
                    if taxon not in taxon_dict[rank]:
                        taxon_dict[rank].append(taxon)
        return taxon_dict
                    
    def calculate_precision(self):
        precision_dict = {'order':[],'family':[], 'genus':[], 'species':[]}
        # true positive / true positive + false positive
        #print(self.correctly_classified)
        for rank in precision_dict.keys():
            correctly_classified_reads_sum = sum ([ c for c in self.correctly_classified[rank]])
            # precision_df = pd.DataFrame(data = {'total_classified_count':total_classified_reads,'unclassified_count':unclassified_reads, 'correctly_classified_reads':correctly_classified_reads})
            # calculating precision with: 
            # true positive = # reads classified to correct taxon (only zymo relevant taxa) 
            # false positive = all classified reads - correctly classified reads + all unclassified reads 
            # final precision calc is correctly_classified_reads / (correctly_classified_reads + (all_classified - correctly_classified))
            precision = correctly_classified_reads_sum / self.total_classified_count
            precision_dict[rank] = precision
        return precision_dict
    
    def calculate_recall(self):
        recall_dict = {'order':[],'family':[], 'genus':[], 'species':[]}
		# recall is tp/(tp+fn)
		# true positive is correctly classified reads 
		# false negative is unclassified reads 
		# final calculation is recall = correctly_classified_reads / (correctly_classified_reads + unclassified_reads)
        for rank in recall_dict.keys():
            correctly_classified_reads_sum = sum ([ c for c in self.correctly_classified[rank]])
            recall = correctly_classified_reads_sum / (self.unclassified + correctly_classified_reads_sum)
            recall_dict[rank] = recall
        return recall_dict
    
    def calculate_f1(self): 
        f1_dict = {'order':[],'family':[], 'genus':[], 'species':[]}
        for rank in f1_dict.keys():
            f1 = (2 * self.precision[rank] * self.recall[rank]) / (self.precision[rank] + self.recall[rank])
            f1_dict[rank] = f1
        return f1_dict
    def calculate_performance_stats(self):
        self.precision = self.calculate_precision()
        self.recall = self.calculate_recall()
        self.f1 = self.calculate_f1()

class read_input():
    def __init__(self, taxon_list, fpath):
        self.taxon_list = taxon_list
        self.fpath = fpath
        #print(self.fpath)

    def readPavian(fpath, s):
        for rank in s.correctly_classified.keys():
            taxon_list = s.taxon_list[rank]
            correctly_classified= [0 for i in range(len(taxon_list))]
            with open(fpath, newline = '\n') as csv_f:
                for row in csv_f:
                    current_taxon = int(row.split()[4])
                    # number of reads at or below the taxa 
                    taxon_count = int(row.split()[1])
                    if current_taxon == 0: 
                        unclassified = taxon_count
                    elif current_taxon == 1:
                        total_classified_count = taxon_count
                    elif current_taxon in taxon_list:
                        correctly_classified[taxon_list.index(current_taxon)] = int(taxon_count)
                s.correctly_classified[rank] = correctly_classified
            s.total_classified_count = total_classified_count
            s.unclassified = unclassified

    def readczID(fpath, upath, s):
        for rank in s.correctly_classified.keys():
            taxon_list = s.taxon_list[rank]
            correctly_classified= [0 for i in range(len(taxon_list))]
            taxon_counts = {}
            nt_values = {}
            with open(fpath, newline = '\n') as csv_f:
                for row in csv_f:
                    if row.startswith(">"):
                        cols = row.split(':')
                        for val in range(len(cols)):
                            if "nt" in val: 
                                r = cols[val].split('_')[0]
                                nt_values[r] = cols[val + 1]
                        current_taxon = nt_values[rank]
                        if current_taxon not in taxon_counts: 
                            taxon_counts[current_taxon] = 1
                        else: 
                            taxon_counts[current_taxon] += 1 
                        # number of reads at or below the taxa 
                unclassified = 0
                unclassified_fasta= open(upath)
                for line in unclassified_fasta:
                    if line.startswith(">"):
                        unclassified += 1
                unclassified_fasta.close()
            for t in taxon_counts.keys:
                if t in taxon_list:
                    correctly_classified[taxon_list.index(t)] = int(taxon_counts[t])
            total_classified_count = sum(taxon_counts.values())
            s.correctly_classified[rank] = correctly_classified
            s.total_classified_count = total_classified_count
            s.unclassified = unclassified

    def readMegan(taxon_list, fpath, name, sample): 
        correctly_classified= [0 for i in range(len(taxon_list))]
        with open(fpath, newline = '\n') as csv_f:
            for row in csv_f:
                current_taxon = int(float(row.split()[0]))
                taxon_count = int(float(row.split()[1]))
                # for megan output unclassified is -1 and -2 is not assigned
                if current_taxon == -1: 
                    unclassified = taxon_count
                elif current_taxon == 1: 
                    total_classified_count = taxon_count
                elif current_taxon == -2: 
                    #total_classified_count = total_classified_count - taxon_count
                    continue
                elif current_taxon in taxon_list: 
                    correctly_classified[taxon_list.index(current_taxon)] = taxon_count
        return(sample(name, unclassified, correctly_classified, total_classified_count))

species_taxid = [1423,5207,1351, 562, 1613, 1639, 287, 4932, 28901, 1280]

zymo_apomp = sample('zymo_apomp', species_taxid)
read_input.readPavian('/Volumes/metagenomics_drive/apomp/publication/zymo_runs/apomp/ERR3152364_GridION_sequencing_EVEN.final.report.tsv', zymo_apomp)
zymo_apomp.calculate_performance_stats()

zymo_mmseq2 = sample('zymo_mmseq2', species_taxid)
read_input.readPavian('/Volumes/metagenomics_drive/apomp/publication/zymo_runs/mmseq2/ERR3152364_GridION_sequencing_EVEN.mmseq.result_report',zymo_mmseq2)
zymo_mmseq2.calculate_performance_stats()

zymo_kraken2 = sample('zymo_kraken2', species_taxid)
read_input.readPavian('/Volumes/metagenomics_drive/apomp/publication/zymo_runs/kraken2/even_out.kraken',zymo_kraken2)
zymo_kraken2.calculate_performance_stats()

zymo_kraken2 = sample('zymo_czid', species_taxid)
read_input.readPavian('/Volumes/metagenomics_drive/apomp/publication/zymo_runs/kraken2/even_out.kraken',zymo_kraken2)
zymo_kraken2.calculate_performance_stats()


#write FPR statistics
fields = ['name','rank','precision','recall','f1']
pipelines = [zymo_apomp, zymo_mmseq2, zymo_kraken2]
#rows = [[n.name for n in pipelines],[p.precision for p in pipelines],[r.recall for r in pipelines],[f.f1 for f in pipelines]]
#rows = [[p.name, p.recall.keys(), p.precision, p.recall, p.f1] for p in pipelines]
rows = []
for p in pipelines:
    if len(rows) == 0:
        r1 = ['pipeline']
        for r in p.recall.keys():
            r1.append('recall_' + str(r))
            r1.append('precision_' + str(r))
            r1.append('f1_' + str(r))
        rows.append(r1)
    temp_row = []
    temp_row.append(p.name)
    print(p.recall)
    for r in p.recall.keys():
        temp_row.append(p.recall[r])
        temp_row.append(p.precision[r])
        temp_row.append(p.f1[r])
    rows.append(temp_row)
filename = "/Volumes/metagenomics_drive/apomp/publication/zymo_runs/analysis/test.csv"
with open(filename, 'w') as csvfile: 
    csvwriter = csv.writer(csvfile) 
    csvwriter.writerows(rows)

# # write per taxa read counts
# fields = [str(t) for t in all_taxids]
# rows = [p.correctly_classified for p in pipelines] 
# filename = "/Volumes/metagenomics_drive/apomp/publication/zymo_runs/analysis/reads_per_taxid.csv"
# with open(filename, 'w') as csvfile: 
#     csvwriter = csv.writer(csvfile) 
#     csvwriter.writerow(fields) 
#     csvwriter.writerows(rows)

