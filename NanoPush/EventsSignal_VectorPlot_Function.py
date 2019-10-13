def SignalVectorsForR(fast5_file):
    """ Produces vectors (tab-separated lists) with signal values (dots), signal average (line)
        & model signals from Nanopolish (rectangle) & my new method = Nanopush (rectangle).
        Signal (y-axis) against contig position (x-axis). Insert file directory with HDF5 file."""

    from AlignQueToRef import Alignment
    call = Alignment(fast5_file, "/Users/kristinaulicna/Documents/Rotation_1/Plotting_Vectors_R/")
    my_list = call.AlignQueToRef()
    contig = call.hit.ctg

    # Dots - signal per each event
    signal_values = []
    for my_line in my_list:
        signal_values.append(float(my_line[3]))

    # Contig positions + How many repetitions per contig position
    ctg_pos_range = list(range(int(my_list[0][0]), int(my_list[-1][0]) + 1))
    how_many = []
    for ctg_pos in ctg_pos_range:
        ctg_pos_counter = 0
        for my_line in my_list:
            if int(my_line[0]) == ctg_pos:
                ctg_pos_counter += 1
        how_many.append(ctg_pos_counter)

    # Contig position + fraction of position depending on the number of repetitions
    ctgpos_values = []
    for repetition, ctg_pos in zip(how_many, ctg_pos_range):
        if repetition != 0:
            for fraction in range(repetition + 1):
                if fraction == 0:
                    continue
                ctgpos_values.append(float(ctg_pos) + (fraction * (1 / (repetition + 1))))

    # K-mers corresponding to each contig:
    ctg_pos_kmer = []
    for repetition, ctg_pos in zip(how_many, ctg_pos_range):
        if repetition == 0:
            ctg_pos_kmer.append(None)
        else:
            for my_line in my_list:
                if ctg_pos == int(my_line[0]):
                    ctg_pos_kmer.append(my_line[8])
                    break

    # Average line for each ctg.position:
    import numpy as np
    mover = 0
    signal_average = []
    for repetition in how_many:
        if mover + repetition <= sum(how_many):
            if repetition == 0:
                signal_average.append(None)
            elif repetition == 1:
                signal_average.append(signal_values[mover])
                mover += 1
            else:
                signal_average.append(np.mean(signal_values[mover:mover+repetition]))
                mover += repetition

    # NANOPOLISH (old) Model rectangle values per each k-mer (mean ± st.dev):
    from Nanopolish_Nanopush_Function import Nanopolish_dictionary
    model_dict_old = Nanopolish_dictionary(kmer_length = '5')
    kmer_values_old = []
    for kmer_old in ctg_pos_kmer:
        for key, value in model_dict_old.items():
            if kmer_old == None:
                kmer_values_old.append([None, None, None])
                break
            elif kmer_old == key:
                kmer_values_old.append(
                    [float(value[0]), float(value[0]) + float(value[1]), float(value[0]) - float(value[1])])

    # NANOPUSH (new) Model rectangle values per each k-mer (mean ± st.dev):
    from Nanopolish_Nanopush_Function import Nanopush_dictionary
    model_dict_new = Nanopush_dictionary()
    kmer_values_new = []
    for kmer_new in ctg_pos_kmer:
        for key, value in model_dict_new.items():
            if kmer_new == None:
                kmer_values_new.append([None, None, None])
                break
            elif kmer_new == key:
                kmer_values_new.append(
                    [float(value[0]), float(value[0]) + float(value[1]), float(value[0]) - float(value[1])])

    """
    # Vectors to return or print:
    print ("Contig:", contig)
    header_1 = "\nI N D I V I D U A L   D O T S (x-axis, y-axis):"
    header_2 = "\nA V E R A G E   O F   D O T S (x-axis, y-axis):"
    header_3 = "\n(OLD) N A N O P O L I S H   M O D E L (x-axis, y-axis):"
    header_4 = "\n(NEW) N A N O P U S H   M O D E L (x-axis, y-axis):"
    print (header_1)
    print (len(ctgpos_values), type(ctgpos_values[0]), ctgpos_values)
    print (len(signal_values), type(signal_values[0]), signal_values)
    print(header_2)
    print(len(ctg_pos_range), type(ctg_pos_range[0]), ctg_pos_range)
    print(len(signal_average), type(signal_average[0]), signal_average)
    print(header_3)
    print(len(ctg_pos_range), type(ctg_pos_range[0]), ctg_pos_range)
    print(len(kmer_values_old), type(kmer_values_old[0]), kmer_values_old)
    print(header_4)
    print(len(ctg_pos_range), type(ctg_pos_range[0]), ctg_pos_range)
    print(len(kmer_values_new), type(kmer_values_new[0]), kmer_values_new)
    """

    """
    # Create a file with these vectors:
    name = fast5_file.split("_")
    name = name[-5] + "_" + name[-4] + "_" + name[-3] + "_" + name[-2]
    file = open("/Users/kristinaulicna/Documents/Rotation_1/Plotting_Vectors_R/" + "/" + "Vectors_" + name + ".txt", "w")

    file.write(header_1+"\n")
    x_val_1 = ''
    y_val_1 = ''
    for item_x, item_y in zip(ctgpos_values, signal_values):
        x_val_1 += (str(item_x) + "\t")
        y_val_1 += (str(item_y) + "\t")
    x_val_1 = x_val_1[:-1]
    y_val_1 = y_val_1[:-1]
    file.write(x_val_1+"\n")
    file.write(y_val_1+"\n")

    file.write(header_2+"\n")
    x_val_2 = ''
    y_val_2 = ''
    for item_x, item_y in zip(ctg_pos_range, signal_average):
        x_val_2 += (str(item_x) + "\t")
        y_val_2 += (str(item_y) + "\t")
    x_val_2 = x_val_2[:-1]
    y_val_2 = y_val_2[:-1]
    file.write(x_val_2+"\n")
    file.write(y_val_2+"\n")

    file.write(header_3+"\n")
    x_val_3 = x_val_2
    y_val_3 = ''
    for lst in kmer_values_old:
        for item in lst:
            y_val_3 += (str(item) + "\t")
    y_val_3 = y_val_3[:-1]
    file.write(x_val_3+"\n")
    file.write(y_val_3+"\n")

    file.write(header_4 + "\n")
    x_val_4 = x_val_2
    y_val_4 = ''
    for lst in kmer_values_new:
        for item in lst:
            y_val_4 += (str(item) + "\t")
    y_val_4 = y_val_4[:-1]
    file.write(x_val_4 + "\n")
    file.write(y_val_4 + "\n")
    file.close()
    """
    return ctgpos_values, signal_values, ctg_pos_range, signal_average, kmer_values_old, kmer_values_new


    # File structure:

# Line[0] = Header #1
# Line[1] = Individual dots (x-axis)
# Line[2] = Individual dots (y-axis)
# Line[3] = b l a n k
# Line[4] = Header #2
# Line[5] = Avg of dots (x-axis)
# Line[6] = Avg of dots (y-axis)
# Line[7] = b l a n k
# Line[8] = Header #3
# Line[9] = Nanopolish Rectangle (x-axis)
# Line[10] = Nanopolish Rectangle (y-axis)
# Line[11] = b l a n k
# Line[12] = Header #4
# Line[13] = Nanopush Rectangle (x-axis)
# Line[14] = Nanopush Rectangle (y-axis)
# Line[15] = b l a n k


    # Call the function:
#ctgpos_values, signal_values, ctg_pos_range, signal_average, kmer_values_old, kmer_values_new = \
#    SignalVectorsForR("/Users/kristinaulicna/Documents/Rotation_1/NanoporeSeq_fast5/CLN_SMD_048987_20180430_FAH86090_MN27963_sequencing_run_BeeRapidNormal_81579_read_6_ch_479_strand.fast5")