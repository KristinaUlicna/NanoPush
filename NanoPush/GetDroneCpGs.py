import time
from WGBSeqClassif_Functions import CreateWGBSeqDict
WGBSeq_dict = CreateWGBSeqDict()

start_time = time.process_time()

unmeth_list = []
meth_list = []
for CpG_line in open("/Users/kristinaulicna/Documents/Rotation_1/Archive/DRONE/NewSelectedFromRob/Sorted_Drone_Sel_Merged.txt", "r"):
    CpG_line = CpG_line.strip().split('\t')
    if CpG_line[6] == 'None' or CpG_line[7] == 'None' or CpG_line[8] == 'None' or CpG_line[9] == 'None' or CpG_line[10] == 'None':
        continue
    else:
        data = [float(number) for number in CpG_line[6:11]]
        CpG_line = CpG_line[0:6] + data
    if CpG_line[0] in WGBSeq_dict[0].keys():
        if int(CpG_line[1]) + 1 in WGBSeq_dict[0][CpG_line[0]]:
            unmeth_list.append(CpG_line)
    if CpG_line[0] in WGBSeq_dict[2].keys():
        if int(CpG_line[1]) + 1 in WGBSeq_dict[2][CpG_line[0]]:
            meth_list.append(CpG_line)


unmeth_file = open("/Users/kristinaulicna/Documents/Rotation_1/Archive/DRONE/NewSelectedFromRob/CpGs_sel_whole_line_Drone_meth_less_0_1.txt", "w")
string = ''
for single_line in unmeth_list:
    for parameter in single_line:
        string += (str(parameter) + "\t")
    string = string[:-1]
    string += "\n"
unmeth_file.write(string)
unmeth_file.close()

meth_file = open("/Users/kristinaulicna/Documents/Rotation_1/Archive/DRONE/NewSelectedFromRob/CpGs_sel_whole_line_Drone_meth_over_0_9.txt", "w")
string = ''
for single_line in meth_list:
    for parameter in single_line:
        string += (str(parameter) + "\t")
    string = string[:-1]
    string += "\n"
meth_file.write(string)
meth_file.close()

print ("File sorting took", time.process_time() - start_time, "seconds, i.e.", (time.process_time() - start_time) / 60, "minutes.")