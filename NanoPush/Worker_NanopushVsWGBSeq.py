#TODO: Create a structure of WGBSeq to be the same as in Nanopolish_reference from Rob:

# 1.) In terminal, sort the file:       sort GSM1360877_worker_brain_10x.txt -k 1,1 -k 2,2n > WGBSeq_Reference.txt

# 2.) Add up the lines covering the same position on + & - strands:

"""
temporary_list = []
WGBSeq_file = open("/Users/kristinaulicna/Documents/Rotation_1/Data/WGBSeq_Reference_Plus_Only.txt", "w")
for counter, line in enumerate(open("/Users/kristinaulicna/Documents/Rotation_1/Data/WGBSeq_Reference_Sorted.txt", "r")):
    if counter > 1000:
        break
    line = line.strip().split('\t')
    line[1] = int(line[1])
    line[4] = int(line[4])
    line[5] = int(line[5])
    if counter == 0:
        temporary_list.append(line)
        continue
    else:
        if line[2] == '+':
            if line[1] + 1 > temporary_list[0][1]:
            temporary_list = []
            temporary_list.append(line)
        if line[2] == '-':
            if line[1] - 1 == temporary_list[0][1]:
                temporary_list[0][4] += line[4]
                temporary_list[0][5] += line[5]
            else:
                temporary_list = []
                temporary_list.append(line)
            string = ''
            string += str(temporary_list[0][0]) + "\t" + str(temporary_list[0][1]) + "\t" + str(
                temporary_list[0][1]) + "\t" + str(temporary_list[0][4] + temporary_list[0][5]) + "\t" + str(
                temporary_list[0][4]) + "\t" + str(
                temporary_list[0][4] / (temporary_list[0][4] + temporary_list[0][5])) + "\n"
            WGBSeq_file.write(string)
WGBSeq_file.close()
print ("Done!")

"""

# 3.) Prepare the file to be the same structure as Rob's Nanopolish:

"""
WGBSeq_file = open("/Users/kristinaulicna/Documents/Rotation_1/Data/WGBSeq_Reference.txt", "w")
for counter, line in enumerate(open("/Users/kristinaulicna/Documents/Rotation_1/Data/WGBSeq_Reference_Sorted.txt", "r")):
    line = line.strip().split('\t')
    line[1] = int(line[1])
    line[4] = int(line[4])
    line[5] = int(line[5])
    WGBSeq_file.write(str(line[0]) + "\t" + str(line[1]) + "\t" + str(line[1]) + "\t" + str(line[4] + line[5]) + "\t" \
                      + str(line[4]) + "\t" + str(line[4] / (line[4] + line[5])) + "\n")
WGBSeq_file.close()

"""

# 4.) Compare the scores given by WGBSeq with my SVM:
"""
from ML_SVM_Worker import TrainSVMonWorker
clf = TrainSVMonWorker()
print (clf)

def WGBSeqVersusSVM(iterable_file, reference_file, reference_file_combo):
    output_file = open(reference_file_combo, "w")
    input_file = open(reference_file, "r")
    CpG = list(input_file.readline().strip().split('\t')) + ['0', '0', '0']
    for line in open(iterable_file, "r"):
        line = line.strip().split('\t')
        data = [float(i) for i in line[6:11]]
        if line[0] == CpG[0]:
            if line[0] == CpG[0] and int(line[1]) == (int(CpG[1])-1):
                CpG[6] = int(CpG[6]) + 1
                if int(clf.predict([data])) == 1:
                    CpG[7] = int(CpG[7]) + 1
            elif line[0] == CpG[0] and int(line[1]) > (int(CpG[1])-1):
                if int(CpG[6]) == 0:
                    CpG[8] = float(CpG[8])
                else:
                    CpG[8] = int(CpG[7]) / int(CpG[6])
                if int(CpG[6]) >= 5:
                    string = ''
                    for parameter in CpG:
                        string += (str(parameter) + "\t")
                    string = string[:-1]
                    string += "\n"
                    output_file.write(string)
                CpG = list(input_file.readline().strip().split('\t')) + ['0', '0', '0']
                if line[0] == CpG[0] and int(line[1]) == (int(CpG[1])-1):
                    CpG[6] = int(CpG[6]) + 1
                    if int(clf.predict([data])) == 1:
                        CpG[7] = int(CpG[7]) + 1
            elif line[0] == CpG[0] and int(line[1]) < (int(CpG[1])-1):
                continue
        else:
            CpG = list(input_file.readline().strip().split('\t')) + ['0', '0', '0']


WGBSeqVersusSVM("/Users/kristinaulicna/Documents/Rotation_1/Data/Sorted_Worker_2_NoNone.txt",\
                    "/Users/kristinaulicna/Documents/Rotation_1/Data/WGBSeq_Reference.txt",\
                    "/Users/kristinaulicna/Documents/Rotation_1/Data/WGBSeq_Reference_SVM.txt")
print ("Done!")
"""

axis_Ref = []
axis_SVM = []
for counter, line in enumerate(open("/Users/kristinaulicna/Documents/Rotation_1/Data/WGBSeq_Reference_SVM.txt", "r")):
    line = line.strip().split('\t')
    axis_Ref.append(float(line[5]))
    axis_SVM.append(float(line[8]))

print (counter)
print (len(axis_Ref), axis_Ref[0:10])
print (len(axis_SVM), axis_SVM[0:10])

import matplotlib.pyplot as plt
plt.scatter(axis_SVM, axis_Ref, alpha=0.2)
plt.xlim(-0.05, 1.05)
plt.ylim(-0.05, 1.05)
plt.xlabel("My Worker-trained SVM")
plt.ylabel("WGBSeq reference")
plt.title("Difference in CpG methylation scoring (coverage > 5)\nwith WGBSeq & my Worker-trained SVM")
#plt.savefig("/Users/kristinaulicna/Documents/Rotation_1/Archive/ScatterNanopolishVsSVM_Drone.png", bbox_inches = 'tight')
plt.show()
