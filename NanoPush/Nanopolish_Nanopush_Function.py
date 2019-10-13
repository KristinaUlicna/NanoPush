def Nanopolish_dictionary(kmer_length = ''):
    model_list = []
    for line in open("/Users/kristinaulicna/Documents/Rotation_1/Nanopolish_" + str(kmer_length) + "mer_model_values.txt", "r"):
        line = line.strip().split('\t')
        if line[0] != 'kmer':
            model_list.append([str(line[0]), float(line[1]), float(line[2])])
    model_dict = {}
    for kmer in model_list:
        model_dict.update({kmer[0]:[kmer[1], kmer[2]]})
    return model_dict


def Nanopush_dictionary():
    Push_dict = {}
    for line in open("/Users/kristinaulicna/Documents/Rotation_1/NanoPush_5mer_Model.txt", "r"):
        line = line.strip().split('\t')
        if str(line[0]) not in Push_dict:
            Push_dict[str(line[0])] = [float(line[1]), float(line[2])]
        Push_dict[str(line[0])] = [float(line[1]), float(line[2])]
    return Push_dict