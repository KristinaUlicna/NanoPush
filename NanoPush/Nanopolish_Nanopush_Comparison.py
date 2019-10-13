#TODO: Visualise the signals for 5-mers using Nanopolish & Nanopush; plot the differences for each 5-mer.

from Nanopolish_Nanopush_Function import Nanopolish_dictionary, Nanopush_dictionary

Nanopolish_5mer_dict = Nanopolish_dictionary(kmer_length='5')
print ("Nanopolish_5mer_dict:", len(Nanopolish_5mer_dict), Nanopolish_5mer_dict)

Push_dict = Nanopush_dictionary()
print ("Nanopush_dict:", len(Push_dict), Push_dict)


    # Plot a lineplot to visualise the differences between Nanopolish & Nanopush signal for a range of 5-mers:

import matplotlib.pyplot as plt
portion_start, portion_end = 105, 140   #range between 100 and 110 works perfectly

# Merge dictionaries:
merged_dict = {}
for key_Pol, value_Pol in Nanopolish_5mer_dict.items():
    for key_Push, value_Push in Push_dict.items():
        if key_Pol == key_Push:
            if key_Pol not in merged_dict:
                merged_dict[key_Pol] = [float(value_Pol[0]),
                                        float(value_Pol[0]) - float(value_Pol[1]),
                                        float(value_Pol[0]) + float(value_Pol[1]),
                                        float(value_Push[0]),
                                        float(value_Push[0]) - float(value_Push[1]),
                                        float(value_Push[0]) + float(value_Push[1])]

print ("\nMerged_Dict:", "\t", len(merged_dict), "\t", merged_dict)


    # Prepare the data:

axes_data = sorted(merged_dict.items())[portion_start:portion_end]
x_axis, y_axis = zip(*axes_data)
x_axis = list(x_axis)
y_axis = list(y_axis)

data = [[] for i in range(6)]
index = 0
while index < 6:
    for mini_list in y_axis:
        data[index].append(mini_list[index])
    index += 1
#print("Datapoints to plot:", data)

print ("X-axis:", len(x_axis), x_axis)
print ("Y-axis:", len(data[0]), data[0])


    # Plot the thing:

plt.figure(figsize = (15, 7))
plt.plot(x_axis, data[0], color='c', marker = 'o', markersize = 10, linewidth = 0.0)
plt.plot(x_axis, data[3], color='coral', marker = 'o', markersize = 10, linewidth = 0.0)
                # I don't want the points to be connected by line.

plt.fill_between(x_axis, data[0], data[1], color='c', alpha=0.3)
plt.fill_between(x_axis, data[0], data[2], color='c', alpha=0.3)
plt.fill_between(x_axis, data[3], data[4], color='coral', alpha=0.3)
plt.fill_between(x_axis, data[3], data[5], color='coral', alpha=0.3)

plt.xticks(fontsize = 12, rotation = 45)
plt.xlabel("K-mer (k = 5, i.e. pentamer)", fontsize = 20, labelpad = 10)
plt.yticks(fontsize = 16)
plt.ylabel("5-mer-specific events signal [pAmp]", fontsize = 20, labelpad = 10)
plt.ylim(60, 135)
plt.grid(b=None, axis='y')

plt.plot([], c='c', marker = 'o', markersize = 10, label='(OLD) Nanopolish signal model (mean ± sd)')
plt.plot([], c='coral', marker = 'o', markersize = 10, label='(NEW) Nanopush signal model (mean ± sd)')
plt.legend(fontsize = 20)

#plt.title("K-mer signal differences between Nanopolish & Nanopush models")
plt.savefig("/Users/kristinaulicna/Documents/Rotation_1/Plots_For_Report/Nanopolish_Nanopush_Signal_Difference.png", bbox_inches = 'tight')
plt.show()
