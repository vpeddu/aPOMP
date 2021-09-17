import sys
import pysam
import logging
from typing import Dict, List, Tuple

LOGGER = logging.getLogger()
LOGGER.setLevel(0)


class Node:
    """ Definition of the class Node """
    def __init__(self):
        self.tax_id = "0"     # Number of the tax id.
        self.parent = "0"     # Number of the parent of this node
        self.children = []    # List of the children of this node
        self.division = None  # Division.
        self.is_tip = True    # Tip = True if it's a terminal node, False if not.
        self.name = ""        # Name of the node: taxa if it's a terminal node, number if not.


def load_ncbi_names(filename: str = "names.dmp") -> Tuple[Dict, Dict]:
    """Load NCBI names definition ("names.dmp")
    Args:
        filename (str): filename of NCBI names
    Returns:
        name_dict, name_dict_reverse
    """

    name_dict = {}  # Initialise dictionary with TAX_ID:NAME
    name_dict_reverse = {}  # Initialise dictionary with NAME:TAX_ID

    LOGGER.warning(f"Load {filename}")
    name_file = open(filename, "r")
    while 1:
        line = name_file.readline()
        if line == "":
            break
        line = line.rstrip()
        line = line.replace("\t", "")
        tab = line.split("|")
        if tab[3] == "scientific name":
            tax_id, name = tab[0], tab[1]  # Assign tax_id and name ...
            name_dict[tax_id] = name  # ... and load them
            name_dict_reverse[name] = str(tax_id)  # ... into dictionaries
    name_file.close()
    return name_dict, name_dict_reverse


def load_ncbi_taxonomy(name_dict, filename: str = "nodes.dmp"):
    """Load taxonomy NCBI file ("nodes.dmp")
    Args:
        filename (str): filename of ncbi taxonomy
        name_dict (dict): name_dict
    Returns:
    """

    # Define taxonomy variable
    # global name_object
    name_object: Dict = {}

    LOGGER.warning(f"Load {filename}")
    taxonomy_file = open(filename, "r")
    while 1:
        line = taxonomy_file.readline()
        if line == "":
            break
        line = line.replace("\t", "")
        tab = line.split("|")

        tax_id = str(tab[0])
        tax_id_parent = str(tab[1])
        division = str(tab[2])

        # Define name of the taxonomy id
        name = "unknown"
        if tax_id in name_dict:
            name = name_dict[tax_id]

        if tax_id not in name_object:
            name_object[tax_id] = Node()
        name_object[tax_id].tax_id = tax_id  # Assign tax_id
        name_object[tax_id].parent = tax_id_parent  # Assign tax_id parent
        name_object[tax_id].name = name  # Assign name
        name_object[tax_id].division = division  # Assign name

        # Add it has children to parents
        children_list = []
        if tax_id_parent in name_object:
            children_list = name_object[
                tax_id_parent
            ].children  # If parent is in the object
        else:
            name_object[tax_id_parent] = Node()
            name_object[tax_id_parent].tax_id = tax_id_parent  # Assign tax_id
        children_list.append(tax_id)  # ... we found its children.
        name_object[
            tax_id_parent
        ].children = children_list  # ... so add them to the parent

        # As the parent node is found, it is not a terminal node then
        name_object[tax_id_parent].is_tip = False

    taxonomy_file.close()
    return name_object

def get_common_ancestor(name_object, node_list: List[str]):
    """
    Function to find common ancestor between two nodes or more
    Args:
        name_object (name_object): taxonomy to use
        node_list (list): list of node
    Returns:
        node (str): node of the common ancestor between nodes
    """

    # global name_object
    list1 = get_genealogy(
        name_object, node_list[0]
    )  # Define the whole genealogy of the first node
    ancestral_list: List[str] = []
    for node in node_list:
        list2 = get_genealogy(
            name_object, node
        )  # Define the whole genealogy of the second node
        ancestral_list = []
        for taxid in list1:
            if taxid in list2:  # Identify common nodes between the two genealogy
                ancestral_list.append(taxid)
        list1 = ancestral_list  # Reassigning ancestral_list to list 1.
    last_common_ancestor = ancestral_list[0]  # Finally, the first node of the ancestral_list is the common ancestor of all nodes.
    return last_common_ancestor  # Return a node

def get_genealogy(name_object, leaf_node: str) -> List[str]:
    """ Trace genealogy from root to leaf """
    ancestors = []  # Initialise the list of all nodes from root to leaf.
    gen_tax_id = leaf_node  # Define leaf
    while 1:
        if gen_tax_id in name_object:
            ancestors.append(gen_tax_id)
            gen_tax_id = name_object[gen_tax_id].parent  # Move up to parents
        else:
            break
        if gen_tax_id == "1":
            # If it is the root, we reached the end.
            # Add it to the list and break the loop
            ancestors.append(gen_tax_id)
            break
    return ancestors  # Return the list

name_dict, name_dict_reverse = load_ncbi_names(filename="names.dmp")  # Load names
ncbi_taxonomy = load_ncbi_taxonomy(filename="nodes.dmp", name_dict=name_dict)
#print(get_common_ancestor(ncbi_taxonomy, ['9606','9601']))

class read():
    def __int__(self):
        self.id = ''
        self.cigar = ''
        self.mapq = ''
        self.seq = ''
        self.taxid = []
        self.seen = false


bamfile = pysam.AlignmentFile(sys.argv[1], "rb")
read_dict = {}
for record in bamfile: 
    #print(record.query_name)
    if record.query_name not in read_dict:
        read_dict[record.query_name] = read()
        read_dict[record.query_name].id = record.query_name
        read_dict[record.query_name].cigar = [record.cigarstring]
        read_dict[record.query_name].mapq = [record.query_qualities]
        read_dict[record.query_name].seq = record.query_sequence
        read_dict[record.query_name].taxid = [record.reference_name]
        read_dict[record.query_name].seen = True
    else: 
        read_dict[record.query_name].cigar.append(record.cigarstring)
        read_dict[record.query_name].mapq.append(record.query_qualities)
        read_dict[record.query_name].taxid.append(record.reference_name)
#print(read_dict['SRR11786979.760008'].cigar)

assignments = {}
for read in read_dict.keys():
    #print(read, read_dict[read].taxid)
    print(read_dict[read].taxid)
    lca = get_common_ancestor(ncbi_taxonomy,read_dict[read].taxid)
    if lca not in assignments:
        assignments[lca] = 1
    else: 
        assignments[lca] += 1

outfilename = sys.argv[2] + '.prekraken.tsv'

with open(outfilename, 'w') as prekraken:
    for taxa in assignments.keys():
        line = str(taxa) + '\t' + str(assignments[taxa])
        prekraken.write("%s\n" % line)




