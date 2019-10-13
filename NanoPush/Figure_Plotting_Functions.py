def EventSignalDotPlot(vector_file_forR):
    vector_list = [[] for i in range(8)]
    for counter, line in enumerate(open(vector_file_forR, "r")):
        line.strip().split('\t')
        if counter == 1:
            for number in line:
                vector_list[0].append(float(number))
        if counter == 2:
            for number in line:
                vector_list[1].append(float(number))
        if counter == 5:
            for number in line:
                vector_list[2].append(float(number))
        if counter == 6:
            for number in line:
                vector_list[3].append(float(number))
        if counter == 9:
            for number in line:
                vector_list[4].append(float(number))
        if counter == 10:
            for number in line:
                vector_list[5].append(float(number))
        if counter == 13:
            for number in line:
                vector_list[6].append(float(number))
        if counter == 14:
            for number in line:
                vector_list[7].append(float(number))
    return vector_list


vector_list = EventSignalDotPlot("/Users/kristinaulicna/Documents/Rotation_1/Plotting_Vectors_R/Vectors_read_6_ch_479.txt")

for index in range(8):
    header = ["\nI N D I V I D U A L   D O T S :", "\nA V E R A G E   O F   D O T S :", "\n(Old) N A N O P O L I S H   M O D E L :", "\n(New) N A N O P U S H   M O D E L :"]
    axis_type = [["x-axis", "y-axis"] * 4]
    if index % 2 == 0:
        print (header[index])
    print (axis_type[index], "\t", len(vector_list[index]), "\t", len(vector_list[index][0]), "\t", vector_list[index])


"""
    # TODO: Plot these in Python!

import matplotlib.pyplot as plt
start = 100
end = 110
index = "[start:end]"

#plt.plot(x_axis, y_axis, "ro")
plt.plot(ctg_pos_range[start:end], signal_average[start:end], "bo")
plt.plot(ctg_pos_range[start:end], kmer_values[start:end], "go")
#plt.xticks(np.arange(min(ctg_pos_range), max(ctg_pos_range)+1, 0.5))
#plt.yticks()
plt.grid(True, color='grey')
plt.show()

"""