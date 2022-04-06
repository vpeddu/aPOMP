import glob

read_files = glob.glob(r'*.unclassified_reads.txt')

s = None
for file in read_files:
	with open(file) as fp:
		if s == None:
			s = set([line.rstrip() for line in fp])
			continue

		new_set = [line.rstrip() for line in fp]
		if len(new_set) != 0:
			s = s.intersection(new_set)

print(len(s))

true_unassigned_file = open("true_unassigned_reads.txt", "w")
print('number of true unassigned reads is ', str(len(s)))
for read in s:
    true_unassigned_file.write(read + "\n")
true_unassigned_file.close()


# import sys
# import math
# import glob

# unassigned_read_files = glob.glob(r'*.unclassified_reads.txt')

# read_id_dict = {}
# for f in unassigned_read_files:
#     with open(f) as fp:
#         Lines = fp.readlines()
#         for line in Lines:
#             if line.strip() not in read_id_dict: 
#                 read_id_dict[line.strip()] = [f]
#             else: 
#                 read_id_dict[line.strip()].append(f)

# true_unassigned = []
# for i in read_id_dict.keys():
#     if len(read_id_dict[i]) == len(unassigned_read_files):
#         true_unassigned.append(i)

# true_unassigned_file = open("true_unassigned_reads.txt", "w")
# print('number of true unassigned reads is ', str(len(true_unassigned)))
# for read in true_unassigned:
#     true_unassigned_file.write(read + "\n")
# true_unassigned_file.close()

#from IPython import embed; embed()

