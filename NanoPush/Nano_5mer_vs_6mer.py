def Model_kmer_Values(Nanopolish_file):
    """Creates a dictionary of model values for each 5-mer.
    Insert the whole name (with whole directory) of the .txt file."""
    model_dict = {}
    model_list = []
    for counter, line in enumerate(open(Nanopolish_file, 'r')):
        if counter != 0:
            parameters = line.strip().split('\t')
            model_list.append(parameters[0:3])
    for k_mer in model_list:
        model_dict.update({k_mer[0]:[k_mer[1], k_mer[2]]})
    return model_list, model_dict

list_5mer, kmer_5_model = Model_kmer_Values("/Users/kristinaulicna/Documents/Rotation_1/Nanopolish_5mer_model_values.txt")
list_6mer, kmer_6_model = Model_kmer_Values("/Users/kristinaulicna/Documents/Rotation_1/Nanopolish_6mer_model_values.txt")

print ("Pentamer Dictionary:", "\t", len(kmer_5_model), "\t", kmer_5_model)
print ("Hexamer Dictionary:", "\t", len(kmer_6_model), "\t", kmer_6_model)


#---------- Working with lists - convert numbers from strings to floats:

for line in list_5mer:
    line[1:3] = [float(number) for number in line[1:3]]
    print (line)

for line in list_6mer:
    line[1:3] = [float(number) for number in line[1:3]]
    print (line)

#----------

divergence_list = []
for line_6 in list_6mer:
    for line_5 in list_5mer:
        if line_5[0] == line_6[0][0:5]:
            divergence_list.append(line_5[1] - line_6[1])
print (len(divergence_list))

import matplotlib.pyplot as plt
plt.boxplot(divergence_list)
plt.title("Distribution of the difference between 5-mer and 6-mer \nNanopolish values for |'C G _ _ _'| position")
plt.xticks()
plt.xlabel("|'C G _ _ _'|")
plt.ylabel("5-mer minus 6-mer")
plt.grid(axis='y')
plt.show()

"""
CGNNN_5_list = []
for key, value in kmer_5_model.items():
    if key[0:2] == "CG":
        CGNNN_5_list.append([key, value])
for line in CGNNN_5_list:
    print (line)

CGNNNN_6_list = []
for key, value in kmer_5_model.items():
    if key[0:2] == "CG":
        for key_2, value_2 in kmer_6_model.items():
"""