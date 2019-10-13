# TODO: Extract unmethylated & methylated sites from the files which Rob will send you.
    # Those are ~200 files each (U & M):
        # The methylated files are expected to have 1 methylated CpG per file.
        # The unmethylated files are expected to have all unmethylated CpGs per file (minus those scored as 10-90%).

#---------- Call the WGBSeq dictionary to which data will be compared:

from WGBSeqClassif_Functions import CreateWGBSeqDict
WGBSeq_dict = CreateWGBSeqDict()
for index in range(3):
    print ("WGBSeq Dict", index, ":", len(WGBSeq_dict[index]), end='\t')
print ()


#---------- Loop through all the files to be analysed (if Rob sends them individually or concatenated, it doesn't matter):

def ExtractCpGsNewModel(folder_directory):
    from os import listdir
    unmeth_list = []
    meth_list = []
    for txt_file in listdir(folder_directory + "/"):
        if txt_file == '.DS_Store':
            #raise Exception("Warning: Remove '.DS_Store' file from the folder in terminal")
            continue
        for CpG_line in open(folder_directory + "/" + txt_file, "r"):
            CpG_line = CpG_line.strip().split('\t')
            if CpG_line[6] == 'None' or CpG_line[7] == 'None' or CpG_line[8] == 'None' or CpG_line[9] == 'None' or CpG_line[10] == 'None':
                continue
            else:
                data = [float(number) for number in CpG_line[6:11]]
            if CpG_line[0] in WGBSeq_dict[0].keys():
                if int(CpG_line[1]) + 1 in WGBSeq_dict[0][CpG_line[0]]:
                    unmeth_list.append(data)
            if CpG_line[0] in WGBSeq_dict[2].keys():
                if int(CpG_line[1]) + 1 in WGBSeq_dict[2][CpG_line[0]]:
                    meth_list.append(data)
    return unmeth_list, meth_list

