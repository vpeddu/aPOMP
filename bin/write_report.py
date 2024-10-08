import sys
import csv
from collections import defaultdict
from ete3 import NCBITaxa
import os

# Initialize NCBITaxa
ncbi = NCBITaxa('taxa.sqlite')

# Read prekraken file
prekraken = sys.argv[1]
taxid_counts = defaultdict(int)

def fix_rank_code(rank):
    rank_to_code = {
        'unclassified': 'U',
        'root': 'R',
        'superkingdom': 'D',
        'kingdom': 'K',
        'phylum': 'P',
        'class': 'C',
        'order': 'O',
        'family': 'F',
        'genus': 'G',
        'species': 'S'
    }
    if rank in rank_to_code:
        return rank_to_code[rank]
    else: 
        return '-'    


def main():
    taxid2node = {}
    with open(prekraken, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:

            taxid_counts[int(row[0])] = int(row[1])
    tree = ncbi.get_topology([key for key in taxid_counts.keys() if key != 0], intermediate_nodes=True)
    print(taxid_counts)
    
    # initialize leafs
    unclassified_count = taxid_counts[0] if 0 in taxid_counts.keys()  else 0
    tree_level_count = taxid_counts[1] if 1 in taxid_counts.keys()  else 0
    tree_root_cum_count = tree_level_count
    for node in tree.traverse("postorder"):
        if node.is_leaf():
            node.add_features(taxon_count= taxid_counts[int(node.name)],
                              cum_sum = taxid_counts[int(node.name)])
            tree_root_cum_count += node.taxon_count

    for node in tree.traverse("postorder"):
        if node.is_leaf():
            #print(node.name)
            print('skipping',node.name)
            continue
        if int(node.name) in taxid_counts.keys():
            children = node.children
            child_sum = 0
            for c in children:
                child_sum += c.cum_sum
            node.add_features(taxon_count= taxid_counts[int(node.name)],
                              cum_sum = taxid_counts[int(node.name)] + child_sum)
        else: 
            child_sum = 0
            children = node.children
            for c in children:
                child_sum += c.cum_sum
            node.add_features(taxon_count= 0,
                              cum_sum = child_sum)
        tree_root_cum_count += node.taxon_count
    for node in tree.traverse("postorder"):
        node.add_feature('percent', 100 * (node.cum_sum/tree_root_cum_count))
        #print(node.name, node.taxon_count, node.cum_sum,node.percent, node.rank, fix_rank_code(node.rank))

    output_file = open(sys.argv[2], 'a')
    output_file.write("  %s" % "{:.2f}".format(float(unclassified_count)/float(tree_root_cum_count)*100))
    output_file.write("\t%i" % unclassified_count)
    output_file.write("\t%i" % unclassified_count)
    output_file.write('\tU\t0\tunclassified\n')
    
    for node in tree.traverse("preorder"):
        space_num = len(node.lineage) - 1
        output_file.write(" %s" % "{:.2f}".format(node.percent))
        output_file.write("\t%i" % node.cum_sum)
        output_file.write("\t%i" % node.taxon_count)
        output_file.write("\t%s" % fix_rank_code(node.rank))
        output_file.write("\t%s" % node.name)
        output_file.write("\t%s" % ("  " * space_num + node.sci_name))
        output_file.write("\n")
    output_file.close()

if __name__ == "__main__":
    main()
