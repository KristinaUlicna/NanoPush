temporary_list = []
WGBSeq_file = open("/Users/kristinaulicna/Documents/Rotation_1/Data/WGBSeq_Reference_Plus_Only.txt", "w")
for counter, line in enumerate(open("/Users/kristinaulicna/Documents/Rotation_1/Data/WGBSeq_Reference_Sorted.txt", "r")):
    line = line.strip().split('\t')
    line[1] = int(line[1])
    line[4] = int(line[4])
    line[5] = int(line[5])
    if counter == 0:
        temporary_list.append(line)
        continue
    else:
        temporary_list.append(line)
        if line[2] == '+':
            if temporary_list[0][2] == '+':
                string = ''
                string += str(temporary_list[1][0]) + "\t" + str(temporary_list[1][1]) + "\t" + str(
                    temporary_list[1][1]) + "\t" + str(temporary_list[1][4] + temporary_list[1][5]) + "\t" + str(
                    temporary_list[1][4]) + "\t" + str(
                    temporary_list[1][4] / (temporary_list[1][4] + temporary_list[1][5])) + "\n"
                WGBSeq_file.write(string)
        if line[2] == '-':
            if line[1] - 1 == temporary_list[0][1]:
                temporary_list[0][4] += line[4]
                temporary_list[0][5] += line[5]
                string = ''
                string += str(temporary_list[0][0]) + "\t" + str(temporary_list[0][1]) + "\t" + str(
                    temporary_list[0][1]) + "\t" + str(temporary_list[0][4] + temporary_list[0][5]) + "\t" + str(
                    temporary_list[0][4]) + "\t" + str(
                    temporary_list[0][4] / (temporary_list[0][4] + temporary_list[0][5])) + "\n"
                WGBSeq_file.write(string)
            else:
                string = ''
                string += str(temporary_list[1][0]) + "\t" + str(temporary_list[1][1]) + "\t" + str(
                    temporary_list[1][1]) + "\t" + str(temporary_list[1][4] + temporary_list[1][5]) + "\t" + str(
                    temporary_list[1][4]) + "\t" + str(
                    temporary_list[1][4] / (temporary_list[1][4] + temporary_list[1][5])) + "\n"
                WGBSeq_file.write(string)
            temporary_list = []
            temporary_list.append(line)

WGBSeq_file.close()
print ("Done!")