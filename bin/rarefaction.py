import os
import sys
import csv
import pandas as pd
import numpy as np
from IPython.display import display

class sample(): 
	def __init__(self, name, taxondict):
		self.name = name
		self.taxondict = taxondict
		self.taxids = self.dict_to_list()[0]
		self.counts = self.dict_to_list()[1]
		self.names = self.dict_to_list()[2]
	def dict_to_list(self):
		tids = []
		cts = []
		names = []
		for i in self.taxondict.keys():
			tids.append(i)
			cts.append(self.taxondict[i])
			names.append(self.name)
		return tids, cts, names
	# def calculate_precision(self):
	# 	# true positive / true positive + false positive
	# 	correctly_classified_reads_sum = sum ([ c for c in self.correctly_classified])
	# 	# precision_df = pd.DataFrame(data = {'total_classified_count':total_classified_reads,'unclassified_count':unclassified_reads, 'correctly_classified_reads':correctly_classified_reads})
	# 	# calculating precision with: 
	# 	# true positive = # reads classified to correct taxon (only zymo relevant taxa) 
	# 	# false positive = all classified reads - correctly classified reads + all unclassified reads 
	# 	# final precision calc is correctly_classified_reads / (correctly_classified_reads + (all_classified - correctly_classified))
	# 	precision = correctly_classified_reads_sum / self.total_classified_count
	# 	return precision

	# def calculate_recall(self):
	# 	# recall is tp/(tp+fn)
	# 	# true positive is correctly classified reads 
	# 	# false negative is unclassified reads 
	# 	# final calculation is recall = correctly_classified_reads / (correctly_classified_reads + unclassified_reads)
	# 	correctly_classified_reads_sum = sum ([ c for c in self.correctly_classified])
	# 	recall = correctly_classified_reads_sum / (self.unclassified + correctly_classified_reads_sum)
	# 	return recall
	# def calculate_f1(self): 
	# 	f1 = (2 * self.precision * self.recall) / (self.precision + self.recall)
	# 	return f1

class read_input():
	def __init__(self, taxon_list, fpath):
		self.taxon_list = taxon_list
		self.fpath = fpath
		print(self.fpath)

	def readPavian(fpath, name): 
		taxondict = {}
		with open(fpath, newline = '\n') as csv_f:
			for row in csv_f:
				current_taxon = int(row.split()[4])
				# number of reads at or below the taxa 
				if row.split()[3] == 'S': 
					taxon_count = int(row.split()[1])
					if current_taxon not in taxondict: 
						taxondict[current_taxon] = taxon_count
		return(sample(name,taxondict))

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

#kraken2 samples
ERR3493775_MinION_sequencing_1 = read_input.readPavian('/Users/vikas/Documents/UCSC/lab/miten_river/ERR3493775_MinION_sequencing_1.final.report.tsv','ERR3493775_MinION_sequencing_1.final.report')
ERR3493774_MinION_sequencing_1 = read_input.readPavian('/Users/vikas/Documents/UCSC/lab/miten_river/ERR3493774_MinION_sequencing_1.final.report.tsv','ERR3493774_MinION_sequencing_1.final.report')
ERR3493787_MinION_sequencing_1 = read_input.readPavian('/Users/vikas/Documents/UCSC/lab/miten_river/ERR3493787_MinION_sequencing_1.final.report.tsv','ERR3493787_MinION_sequencing_1.final.report')
ERR3493789_MinION_sequencing_1 = read_input.readPavian('/Users/vikas/Documents/UCSC/lab/miten_river/ERR3493789_MinION_sequencing_1.final.report.tsv','ERR3493789_MinION_sequencing_1.final.report')
ERR3493793_MinION_sequencing_1 = read_input.readPavian('/Users/vikas/Documents/UCSC/lab/miten_river/ERR3493793_MinION_sequencing_1.final.report.tsv','ERR3493793_MinION_sequencing_1.final.report')
ERR3493777_MinION_sequencing_1 = read_input.readPavian('/Users/vikas/Documents/UCSC/lab/miten_river/ERR3493777_MinION_sequencing_1.final.report.tsv','ERR3493777_MinION_sequencing_1.final.report')

fields = ['name','taxid','count']
pipelines = [ERR3493775_MinION_sequencing_1, ERR3493774_MinION_sequencing_1, ERR3493787_MinION_sequencing_1,ERR3493789_MinION_sequencing_1,ERR3493793_MinION_sequencing_1,ERR3493777_MinION_sequencing_1]
rows = [[p.name, p.taxids, p.counts] for p in pipelines]
filename = "kraken2_raraefaction_stats.csv"
with open(filename, 'w') as csvfile: 
    csvwriter = csv.writer(csvfile) 
    #csvwriter.writerow(fields) 
    #csvwriter.writerows(rows)
    csvwriter.writerows(list(zip(*rows)))

#apomp samples 
fields = ['name','taxid','count']
pipelines = [ERR3493775_MinION_sequencing_1, ERR3493774_MinION_sequencing_1, ERR3493787_MinION_sequencing_1,ERR3493789_MinION_sequencing_1,ERR3493793_MinION_sequencing_1,ERR3493777_MinION_sequencing_1]
rows = [[p.name, p.taxids, p.counts] for p in pipelines]
filename = "kraken2_raraefaction_stats.csv"
with open(filename, 'w') as csvfile: 
    csvwriter = csv.writer(csvfile) 
    #csvwriter.writerow(fields) 
    #csvwriter.writerows(rows)
    csvwriter.writerows(list(zip(*rows)))
