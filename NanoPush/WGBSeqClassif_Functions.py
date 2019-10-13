def CreateWGBSeqDict():
    """Returns single dictionary with 3 subdicts with contig positions according to methylation percentage."""
    WGBSeq_dict = [{} for i in range(3)]
    counter = [0 for i in range(3)]
    for line in open("/Users/kristinaulicna/Documents/Rotation_1/GSM1360877_worker_brain_10x.txt", "r"):
        CpG = line.strip().split('\t')
        #In the WGBSeq file, 5th column (CpG[4]) = methylated and in the 6th column (CpG[5]) = unmethylated. Therefore:
        CpG[4] = (int(CpG[4]) * 100) / (int(CpG[4]) + int(CpG[5]))
        CpG = CpG[:-1]
        if float(CpG[4]) <= 10.0:
            index = 0
        if float(CpG[4]) > 10.0 and float(CpG[4]) < 90.0:
            index = 1
        if float(CpG[4]) >= 90.0:
            index = 2
        counter[index] += 1
        if CpG[0] not in WGBSeq_dict[index]:
            WGBSeq_dict[index][CpG[0]] = [int(CpG[1])]
        else:
            WGBSeq_dict[index][CpG[0]].append(int(CpG[1]))
    print ("\nCpGs in 3 sublayers ([0] = <0.1, [1] = 0.1-0.9, [2] = >0.9):", counter[0], counter[1], counter[2], "\tTotal CpGs:", counter[0] + counter[1] + counter[2])
    return WGBSeq_dict


def ClassifyCpGsWGBSeq(noNone_file, save_to_directory, how_many_unmeth):
    WGBSeq_dict = CreateWGBSeqDict()
    names = ["less_0_1", "bet10_90", "over_0_9"]
    files = [open(save_to_directory + "/CpG_all_Drone_meth_" + name + ".txt", "w") for name in names]
    counters = [0 for i in range(3)]
    excluded_CpG_counter = 0
    for line in open(noNone_file, "r"):
        CpG = line.strip().split('\t')
        index = 3
        #Add one to contig position (int(CpG[1])+1) bc my class counts from 0 while WGBSeq file counts from 1:
        if counters[0] < how_many_unmeth:
            if CpG[0] in WGBSeq_dict[0].keys():
                if (int(CpG[1]) + 1) in WGBSeq_dict[0][CpG[0]]:
                    index = 0
        if CpG[0] in WGBSeq_dict[1].keys():
            if (int(CpG[1]) + 1) in WGBSeq_dict[1][CpG[0]]:
                index = 1
        if CpG[0] in WGBSeq_dict[2].keys():
            if (int(CpG[1]) + 1) in WGBSeq_dict[2][CpG[0]]:
                index = 2
        if index != 0 and index != 1 and index != 2:
            excluded_CpG_counter += 1
            continue
        counters[index] += 1
        if counters[0] % 500 == 0:
            print ("Counters:", names, counters)
        string = ''
        for item in CpG:
            string += item + "\t"
        string = string[:-1]
        string += "\n"
        files[index].write(string)
    [file.close() for file in files]
    print ("File lines count:\tless_0_1 =", counters[0], "\tbet10_90 =", counters[1], "\tover_0_9 =", counters[2])
    print ("Total CpGs:", counters[0] + counters[1] + counters[2] + excluded_CpG_counter, "\tExcluded CpGs:", excluded_CpG_counter)


#ClassifyCpGsWGBSeq("/Users/kristinaulicna/Documents/Rotation_1/Data/Sorted_Drone_NoNone.txt", "/Users/kristinaulicna/Documents/Rotation_1/Data/", 500)