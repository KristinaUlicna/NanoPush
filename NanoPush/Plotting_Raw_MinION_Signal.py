from EventsSignal_VectorPlot_Function import SignalVectorsForR
ctgpos_values, signal_values, ctg_pos_range, signal_average, kmer_values_old, kmer_values_new = \
    SignalVectorsForR("/Users/kristinaulicna/Documents/Rotation_1/NanoporeSeq_fast5/CLN_SMD_048987_20180430_FAH86090_MN27963_sequencing_run_BeeRapidNormal_81579_read_6_ch_479_strand.fast5")

# 1 = Individual points.
# 2 = Average of points.
# 3 = My Nanopush range.
# 4 = Nanopolish range.
# ______________________


    # Section of contig to be visualised:

start = 2 * 120
end = 2 * 130


    # 2

print ("\nA V E R A G E   O F   D O T S (x-axis, y-axis):")

x_2 = ctg_pos_range
x_2 = [[item+0.001, item+0.999] for item in x_2]
x_2 = [j for i in x_2 for j in i]
x_2 = x_2[start:end]
print (len(x_2), x_2)

y_2 = signal_average
y_2 = [item for item in y_2 for i in range(2)]
y_2 = y_2[start:end]
print (len(y_2), y_2)


    # 3

print ("\n(NEW) N A N O P U S H   M O D E L (x-axis, y-axis):")

x_3 = x_2
print (len(x_3), x_3)

y_3 = [[] for i in range(3)]
index = 0
while index < 3:
    for mini_list in kmer_values_new:
        y_3[index].append(mini_list[index])
    index += 1
y_3 = [[item for item in y_3[0] for i in range(2)], [item for item in y_3[1] for i in range(2)], [item for item in y_3[2] for i in range(2)]]
y_3[0] = y_3[0][start:end]
y_3[1] = y_3[1][start:end]
y_3[2] = y_3[2][start:end]
print (len(y_3[0]), len(y_3[1]), len(y_3[2]), y_3)


    # 4

print ("\n(OLD) N A N O P O L I S H   M O D E L (x-axis, y-axis):")

x_4 = x_2
print (len(x_4), x_4)

y_4 = [[] for i in range(3)]
index = 0
while index < 3:
    for mini_list in kmer_values_old:
        y_4[index].append(mini_list[index])
    index += 1
y_4 = [[item for item in y_4[0] for i in range(2)], [item for item in y_4[1] for i in range(2)], [item for item in y_4[2] for i in range(2)]]
y_4[0] = y_4[0][start:end]
y_4[1] = y_4[1][start:end]
y_4[2] = y_4[2][start:end]
print (len(y_4[0]), len(y_4[1]), len(y_4[2]), y_4)


    # 1

print ("\nI N D I V I D U A L   D O T S (x-axis, y-axis):")

min_ctg_pos = int(x_2[0])
print ("Min:", min_ctg_pos)
max_ctg_pos = int(x_2[len(x_2) - 1]) + 1
print ("Max:", max_ctg_pos)

min_index = 0
max_index = 0
for number in ctgpos_values:
    if number >= float(min_ctg_pos):
        min_index = ctgpos_values.index(number)
        break
print ("MIN INDEX:", min_index)

for number in ctgpos_values:
    if number <= float(max_ctg_pos):
        max_index = ctgpos_values.index(number)
print ("MAX INDEX:", max_index)


x_1 = ctgpos_values[min_index:max_index]
print (len(x_1), x_1)

y_1 = signal_values[min_index:max_index]
print (len(y_1), y_1)


    # Ticks:
ticks = []
for position in x_2:
    ticks.append(int(position))
ticks = list(set(ticks))
print ("TICKS:", len(ticks), ticks)

ticks_label = [str(item) for item in ticks]
ticks_posit = [item + 0.5 for item in ticks]


    # Plot the thing:

import matplotlib.pyplot as plt

fig = plt.figure(figsize=(15,7))
ax = plt.subplot(111)

plt.plot(x_3, y_3[0], linewidth = 3.0, color = 'orange')
plt.fill_between(x_3, y_3[0], y_3[1], color = 'orange', alpha = 0.2)
plt.fill_between(x_3, y_3[0], y_3[2], color = 'orange', alpha = 0.2)

plt.plot(x_4, y_4[0], linewidth = 3.0, color = 'green')
plt.fill_between(x_4, y_4[0], y_4[1], color = 'green', alpha = 0.2)
plt.fill_between(x_4, y_4[0], y_4[2], color = 'green', alpha = 0.2)

plt.plot(x_2, y_2, '-', linewidth = 7.0, color = 'dodgerblue')
plt.plot(x_1, y_1, 'ro', markersize = 10, alpha = 0.6)

plt.xticks(ticks_posit, ticks_label, rotation = 45, fontsize = 16)
plt.yticks(fontsize = 16)
plt.ylim(55, 118)

plt.plot([], c='red', marker = 'o', markersize = 10, linestyle = '', alpha = 0.6, label='Individual events \nper contig position')
plt.plot([], c='dodgerblue', linewidth = 7.0, label='Signal average \nper contig position')
plt.plot([], c='green', linewidth = 3.0, label='(OLD) NanoPolish \nk-mer model value \n(mean ± st.dev)')
plt.plot([], c='orange', linewidth = 3.0, label='(NEW) NanoPush \nk-mer model value \n(mean ± st.dev)')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 20)

#plt.title("Raw signal from MinION & segmentation to events")
plt.grid(b=None, which='major', axis='y')
plt.xlabel("Position on contig NC_007075.3 [bp]", fontsize = 20, labelpad = 10)
plt.ylabel("Event signal [pAmp]", fontsize = 20)
plt.savefig("/Users/kristinaulicna/Documents/Rotation_1/Plots_For_Report/Events_Segmentation_Zoom.png", bbox_inches = 'tight')

plt.show()
plt.close()

#TODO: (normalised 5 datapoints average) => to the legend.