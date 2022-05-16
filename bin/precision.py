import os
import sys
import csv
import pandas as pd
import numpy as np
from IPython.display import display

all_taxids = [1423,5207,1351, 562, 1613, 1639, 287, 4932, 28901, 1280]

class sample(): 
	def __init__(self, name, unclassified, correctly_classified, total_classified_count):
		self.name = name
		self.unclassified = unclassified
		self.correctly_classified = correctly_classified
		self.total_classified_count = total_classified_count
		self.precision = self.calculate_precision()
		self.recall = self.calculate_recall()
		self.f1 = self.calculate_f1()

	def calculate_precision(self):
		# true positive / true positive + false positive
		print(self.correctly_classified)
		correctly_classified_reads_sum = sum ([ c for c in self.correctly_classified])
		# precision_df = pd.DataFrame(data = {'total_classified_count':total_classified_reads,'unclassified_count':unclassified_reads, 'correctly_classified_reads':correctly_classified_reads})
		# calculating precision with: 
		# true positive = # reads classified to correct taxon (only zymo relevant taxa) 
		# false positive = all classified reads - correctly classified reads + all unclassified reads 
		# final precision calc is correctly_classified_reads / (correctly_classified_reads + (all_classified - correctly_classified))
		precision = correctly_classified_reads_sum / self.total_classified_count
		return precision

	def calculate_recall(self):
		# recall is tp/(tp+fn)
		# true positive is correctly classified reads 
		# false negative is unclassified reads 
		# final calculation is recall = correctly_classified_reads / (correctly_classified_reads + unclassified_reads)
		correctly_classified_reads_sum = sum ([ c for c in self.correctly_classified])
		recall = correctly_classified_reads_sum / (self.unclassified + correctly_classified_reads_sum)
		return recall
	def calculate_f1(self): 
		f1 = (2 * self.precision * self.recall) / (self.precision + self.recall)
		return f1

class read_input():
	def __init__(self, taxon_list, fpath):
		self.taxon_list = taxon_list
		self.fpath = fpath
		print(self.fpath)

	def readPavian(taxon_list, fpath, name):
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
					correctly_classified[taxon_list.index(current_taxon)] = taxon_count
		return(sample(name, unclassified, correctly_classified, total_classified_count))

	def readMegan(taxon_list, fpath, name): 
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
					total_classified_count = total_classified_count - taxon_count
				elif current_taxon in taxon_list: 
					correctly_classified[taxon_list.index(current_taxon)] = taxon_count
		return(sample(name, unclassified, correctly_classified, total_classified_count))




# desktop
#zymo_guppy_fast_apomp = readpavian(all_taxids,'/media/vikas/thiccy1/data/miten_metagenomics/results_from_mustard/Zymo-GridION-EVEN-BB-SN_Guppy_6.0.1_fast.final.report.tsv')
#zymo_guppy_hac_apomp = readpavian(all_taxids,'/media/vikas/thiccy1/data/miten_metagenomics/results_from_mustard/Zymo-GridION-EVEN-BB-SN_Guppy_6.0.1_hac.final.report.tsv')

#local 
zymo_guppy_sup_apomp_previous = read_input.readPavian(all_taxids,'/Users/vikas/Documents/UCSC/lab/ugh/original_zymo_sup_apomp.tsv','apomp_previous')
zymo_guppy_sup_apomp_plasmid_refactor = read_input.readPavian(all_taxids,'/Users/vikas/Downloads/Zymo-GridION-EVEN-BB-SN_Guppy_6.0.1_sup.final_ntasm20_full.report.tsv','apomp_plasmid_refactor')
zymo_guppy_sup_kraken2 = read_input.readPavian(all_taxids,'/Users/vikas/Downloads/Zymo-GridION-EVEN-BB-SN_Guppy_6.0.1_sup.fastq.gz.kraken2.report','kraken2')
zymo_guppy_sup_megan = read_input.readMegan(all_taxids,'/Users/vikas/Downloads/zymo_sup-ex.txt','diamond_meganlr')

print(zymo_guppy_sup_apomp_previous.total_classified_count)
print(zymo_guppy_sup_apomp_plasmid_refactor.total_classified_count)

fields = ['name','precision','recall','f1']
pipelines = [zymo_guppy_sup_apomp_previous,zymo_guppy_sup_apomp_plasmid_refactor , zymo_guppy_sup_kraken2, zymo_guppy_sup_megan]
#rows = [[n.name for n in pipelines],[p.precision for p in pipelines],[r.recall for r in pipelines],[f.f1 for f in pipelines]]
rows = [[p.name, p.precision, p.recall, p.f1] for p in pipelines]
filename = "zymo_reporting_stats.csv"
with open(filename, 'w') as csvfile: 
    csvwriter = csv.writer(csvfile) 
    csvwriter.writerow(fields) 
    csvwriter.writerows(rows)
