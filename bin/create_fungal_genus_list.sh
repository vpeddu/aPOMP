#!/bin/bash/
# 4751 is the ncbi taxid for fungi
# ete3 leaves out candida and some other big ones for whatever reason 
taxonkit list --ids 4751 | taxonkit filter -E Genus > fungi_genera_list.txt