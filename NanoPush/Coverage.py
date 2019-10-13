# Calculating the coverage of each CpG dinucleotide in the cat-sorted file:
        # TODO: SIZE OF FILE WITH COVERAGE EXTREMELY LOW! (12MB vs ~800MB previously)


# Import entire concatenated and sorted list of CpG dinucleotides from the Maxifile:

maxifile = "/Users/kristinaulicna/Documents/Rotation_1/CpGsOutputMaxifile.txt"
with open(maxifile, 'r') as tsv:
    dinucl_list = []
    for line in tsv:
        parameters = line.strip().split('\t')
        dinucl_list.append(parameters)

print ("Total CpGs in Maxifile: " + str(len(dinucl_list)) + "\t" + "File: " + maxifile.split("/")[5])
#for line in dinucl_list:
#    print (line)


# Filter only those CpG dinucleotides which are covered more than once in the Maxifile:

coverage_list = []
line_order = 0
while line_order < (len(dinucl_list) - 1):
    if dinucl_list[line_order][0] == dinucl_list[line_order+1][0] and dinucl_list[line_order][1] == dinucl_list[line_order+1][1] and dinucl_list[line_order][3] != dinucl_list[line_order+1][3]:
        coverage_list.append(dinucl_list[line_order])
        coverage_list.append(dinucl_list[line_order+1])
    line_order += 1

coverage_list = [coverage_list[i] for i in range(len(coverage_list)) if i == 0 or coverage_list[i] != coverage_list[i - 1]]     # = like a set() function to deduplicate list of lists
print ("\nCoverage List - number of lines printed overall:" + "\t" + str(len(coverage_list)))
for dinucl in coverage_list:
    print (dinucl)


# Calculate the number of unique CpG dinucleotides covered more than once:

line_order = 0
unique_dinucl = 1
        #TODO: unique_dinucl = 1 ??? Is this correct?
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
print ("\nCoverage dictionary with coverage:", "\t", coverage_dict)
print ("\nSummary of coverages per CpG:", "\t", len(coverage_summary), "\t", coverage_summary)


"""
# Access individual layers of the nested dictionary:

for key_1, value_1 in coverage_dict.items():
    print ("Iterating over dictionary:", "\t", key_1, value_1)
    print ("Iterating over k-mers only:", "\t", value_1.keys())
    print ("Iterating over signals only:", "\t", value_1.values())
    break
print ()

"""