# Function to create data sets for training, validation and testing of SVM:

def FilterNone(input_file, output_file):
    noNone_file = open(output_file, "w")
    for CpG_line in open(input_file, "r"):
        string = ''
        CpG_data = CpG_line.strip().split("\t")
        if len(CpG_data) != 11:
            continue
        if CpG_data[6] != 'None' and \
                CpG_data[7] != 'None' and \
                CpG_data[8] != 'None' and \
                CpG_data[9] != 'None' and \
                CpG_data[10] != 'None':
            for parameter in CpG_data:
                string += parameter + "\t"
            string = string[:-1]
            string += "\n"
            noNone_file.write(string)
    noNone_file.close()


def CreateDataSet(input_file, start_line, finish_line):
    CpG_list = []
    for counter, CpG in enumerate(open(input_file, "r")):
        if counter < start_line - 1:
            continue
        if counter > finish_line - 1:
            break
        CpG = CpG.strip().split('\t')
        if len(CpG) == 5:
            CpG = [float(i) for i in CpG]
        elif len(CpG) == 11:
            CpG = [float(i) for i in CpG[6:11]]
        else:
            raise Exception("Warning: Input file has more/less columns than expected")
        CpG_list.append(CpG)
    return CpG_list


def RandomiseLists(unmeth_list, meth_list):
    import random
    merged_list = unmeth_list + meth_list
    meth_score = [0 for i in range(len(unmeth_list))] + [1 for i in range(len(meth_list))]
    scored_CpGs = list(zip(merged_list, meth_score))
    random.shuffle(scored_CpGs)
    merged_list, meth_score = zip(*scored_CpGs)
    merged_list = list(merged_list)
    meth_score = list(meth_score)
    return merged_list, meth_score


def BoxplotPlottingData(input_unmeth_list, input_meth_list):
    from itertools import zip_longest
    unmeth_pos = [[] for i in range(5)]
    meth_pos = [[] for i in range(5)]
    for line_counter, (CpG_un, CpG_me) in enumerate(zip_longest(input_unmeth_list, input_meth_list, fillvalue=None)):
        CpG_counter = 0
        if line_counter <= len(input_unmeth_list):
            while CpG_counter < len(CpG_un):
                unmeth_pos[CpG_counter].append(CpG_un[CpG_counter])
                CpG_counter += 1
        CpG_counter = 0
        if line_counter <= len(input_meth_list) and CpG_me != None:
            while CpG_counter < len(CpG_me):
                meth_pos[CpG_counter].append(CpG_me[CpG_counter])
                CpG_counter += 1
    return unmeth_pos, meth_pos


def NormaliseData(input_list, normalisation_factor):
    for single_list in input_list:
        for order, data_point in enumerate(single_list):
            single_list[order] = data_point + float(normalisation_factor)
    return input_list
