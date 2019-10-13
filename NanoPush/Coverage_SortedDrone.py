# TODO: Calculate the coverage of each CpG dinucleotide in the cat-sorted_dronefile:

# Import entire concatenated and sorted list of CpG dinucleotides from the cat-sorted_drone file:

import time
start_time = time.process_time()

dinucl_list = []
for counter, line in enumerate(open("/Users/kristinaulicna/Documents/Rotation_1/Sorted_Drone_CpGs.txt", "r")):
    dinucl_list.append(line.strip().split('\t'))
print ("\nTotal CpGs in Sorted_Drone_CpGs.txt: " + str(len(dinucl_list)))

# Filter only those CpG dinucleotides which are covered more than once in the Maxifile:

coverage_list = []
line_order = 0
while line_order < (len(dinucl_list) - 1):
    if dinucl_list[line_order][0] == dinucl_list[line_order+1][0] and dinucl_list[line_order][1] == dinucl_list[line_order+1][1] and dinucl_list[line_order][3] != dinucl_list[line_order+1][3]:
        coverage_list.append(dinucl_list[line_order])
        coverage_list.append(dinucl_list[line_order+1])
    line_order += 1
coverage_list = [coverage_list[i] for i in range(len(coverage_list)) if i == 0 or coverage_list[i] != coverage_list[i - 1]]     # = like a set() function to deduplicate list of lists
#print ("\nCoverage List - number of lines printed overall:" + "\t" + str(len(coverage_list)))
#for dinucl in coverage_list:
#    print (dinucl)


# Calculate the number of unique CpG dinucleotides covered more than once:

line_order = 0
unique_dinucl = 1
while line_order < (len(coverage_list) - 1):
    if coverage_list[line_order][0] != coverage_list[line_order + 1][0] or coverage_list[line_order][1] != coverage_list[line_order + 1][1]:
        unique_dinucl += 1
    line_order += 1
print("\nNumber of unique CpG dinucleotides covered more than once:" + "\t" + str(unique_dinucl))


# Make a layered dictionary with a contig as primary key, contig_position as secondary key, and readID as value:

coverage_dict = {}
for dinucl in coverage_list:
    if str(dinucl[0]) not in coverage_dict:
        coverage_dict[str(dinucl[0])] = {}
    if int(dinucl[1]) not in coverage_dict[str(dinucl[0])]:
        coverage_dict[str(dinucl[0])][int(dinucl[1])] = []
    coverage_dict[str(dinucl[0])][int(dinucl[1])].append(str(dinucl[3]))


# Add the coverage number to the end of the list of readIDs which cover this CpG:

coverage_summary = []
for keys, values in coverage_dict.items():
    for key, value in values.items():
        coverage_dict[keys][key].append(len(coverage_dict[keys][key]))
        coverage_summary.append(int(len(coverage_dict[keys][key]) - 1))

coverage_summary = sorted(coverage_summary, reverse=True)
print ("\nMaximum coverage achieved:", coverage_summary[0])
print ("\nMean coverage per CpG covered >2:", sum(coverage_summary)/len(coverage_summary))

print ("\nT I M E :")
print ("Processing took", time.process_time() - start_time, "seconds, i.e.", (time.process_time() - start_time) / 60, "minutes.")