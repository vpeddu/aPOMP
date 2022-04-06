import sys
import csv
import pandas as pd
from IPython.display import display


# keys are taxid, value is list [taxid, name]
true_mapped_read_counts	=	{}
with open('/media/vikas/thiccy1/data/apomp_testing/mapped_read_counts.txt', newline='\n') as csv_f:
	for row in csv_f:
		true_mapped_read_counts[int(row.split()[1])] = [row.split()[1],row.split()[0]]
#print(true_mapped_read_counts)

all_taxids = [int(i) for i in true_mapped_read_counts.keys()]

def readpavian(taxon_list, fpath):
	out_dict={'unclassified_count':'','classified_count':'','taxon_counts':[0 for i in range(len(taxon_list))]}
	c = 0
	with open(fpath, newline = '\n') as csv_f:
		for row in csv_f:
			current_taxon = int(row.split()[4])
			if current_taxon == 0: 
				out_dict['unclassified_count'] = row.split()[1]
			elif current_taxon == 1:
				out_dict['classified_count'] = row.split()[1]
			elif current_taxon in taxon_list:
				out_dict['taxon_counts'][taxon_list.index(current_taxon)] = row.split()[1]
	#print(taxon_list)
	return(out_dict)

zymo_guppy_fast_apomp = readpavian(all_taxids,'/media/vikas/thiccy1/data/miten_metagenomics/results_from_mustard/Zymo-GridION-EVEN-BB-SN_Guppy_6.0.1_fast.final.report.tsv')
zymo_guppy_hac_apomp = readpavian(all_taxids,'/media/vikas/thiccy1/data/miten_metagenomics/results_from_mustard/Zymo-GridION-EVEN-BB-SN_Guppy_6.0.1_hac.final.report.tsv')
zymo_guppy_sup_apomp = readpavian(all_taxids,'/media/vikas/thiccy1/data/miten_metagenomics/results_from_mustard/Zymo-GridION-EVEN-BB-SN_Guppy_6.0.1_sup.final.report.tsv')

def precision(dict_list):
	classified_reads = [int(ccount['classified_count']) for ccount in dict_list]
	unclassified_reads = [int(ucount['unclassified_count']) for ucount in dict_list]
	recall_df = pd.DataFrame(data = {'classified_count':classified_reads,'unclassified_count':unclassified_reads})
	recall_df['recall'] = recall_df['classified_count'] / (recall_df['classified_count'] + recall_df['unclassified_count'])
	display(recall_df)

def recall(dict_list):
	# recall is tp/(tp+fn)
	classified_reads = [int(ccount['classified_count']) for ccount in dict_list]
	unclassified_reads = [int(ucount['unclassified_count']) for ucount in dict_list]
	recall_df = pd.DataFrame(data = {'classified_count':classified_reads,'unclassified_count':unclassified_reads})
	recall_df['recall'] = recall_df['classified_count'] / (recall_df['classified_count'] + recall_df['unclassified_count'])
	display(recall_df)
	


precision([zymo_guppy_hac_apomp,zymo_guppy_fast_apomp, zymo_guppy_sup_apomp]) 