# General directory for HDF5 files:
"/Users/kristinaulicna/Documents/Rotation_1_Rob_Lowe/Nanopore sequences(fast5)/..."


# UNIVERSAL FUNCTION FOR THE ENTIRE ANALYSIS OF ALL (+ve or -ve) HDF5 FILES:

def AlignQueToRef(HDF5_file):
    """Insert file name into "" with a full directory (from /Users/...)
    1.) Extracts the FastQ sequence (string) from an HDF5 file type.
    2.) Creates the alignment of your query sequence to the reference.
    3.) Checks for strand orientation. Leaves positive strand untouched.
    4.) Converts a negative strand into positive by creating complement.
    5.) Re-aligns the optimised strand to the reference & prints both sequences.
    6.) Maps the enumeration of respective mismatches, insertions and deletions
    from the Cigar sequence + prints the analysis:
        - Breaks down of aligned Ref & Que sequences accoring to the Cigar Parameters.
        - Produces a continuous Ref & Que sequences (with respective in/del gaps).
        - Summarizes the total number of matches/mismatches, insertions and deletions.
        - Prints out each individual base in Ref and Que on contig position (no gaps in Ref seq).
        - Checks & prints out mismatches in the aligned sequences.
    RETURNS hit, reference, query_seq, reference_with_gaps, query_with_gaps."""
    print ("HDF5 file name:", "\t", HDF5_file)

    #Extraction of FastQ sequence (string) from an HDF5 file type:
    import h5py
    file = h5py.File(HDF5_file, "r")
    import numpy as np
    readSeq = np.array(file["Analyses"]["Basecall_1D_000"]["BaseCalled_template"]["Fastq"])
    readSeq = str(readSeq)
    readSeq = readSeq.split("\\n")
    readSeq = readSeq[1]
    print ("ReadSeq:", "\t", len(readSeq), "\t", readSeq, "\t", type(readSeq))

    #Alignment of the readSeq sequence to the reference genome (honeybee - Apis mellifera):
    import mappy as mp
    alignment = mp.Aligner("/Users/kristinaulicna/Documents/Rotation_1_Rob_Lowe/GCF_000002195.4_Amel_4.5_genomic.fa.gz",
                           preset="map-ont")
    for hit in alignment.map(readSeq):
        print ("Hit parameters:", "\t", hit)

    #Convertion of the strand (readSeq) into positive orientation:
    if hit.strand == -1:
        query = readSeq[::-1]
        query.upper()
        query = str(query) \
            .replace('T', '%temp1%').replace('A', 'T').replace('%temp1%', 'A') \
            .replace('C', '%temp2%').replace('G', 'C').replace('%temp2%', 'G')
        print ("Strand in negative orientation. Converted seq:", query)
    else:
        query = readSeq
        query.upper()
        print ("Strand in positive orientation. Unchanged seq:", query)
    print ()

    #Change hit.q_st and hit.q_en values if strand was negative:
    if hit.strand == -1:
        reference = alignment.seq(hit.ctg, hit.r_st, hit.r_en)
        query_seq = query[len(query)-hit.q_en:len(query)-hit.q_st]
    else:
        reference = alignment.seq(hit.ctg, hit.r_st, hit.r_en)
        query_seq = query[hit.q_st:hit.q_en]
    print ("Reference:", "\t", len(reference), "\t", reference)
    print ("Query seq:", "\t", len(query_seq), "\t", query_seq)
    print ()

    #Breaking down of aligned Ref & Que sequences accoring to the Cigar Parameters:
    print ("Alignment Analysis - Chunks with matches/mismatches, insertions and deletions:")
    position_ref = 0
    position_que = 0
    for parameter in hit.cigar:
        if parameter[1] == 0:
            print("Ref:", position_ref + 1, ":", position_ref + parameter[0], "\t",
                  reference[position_ref:position_ref + parameter[0]])
            print("Que:", position_que + 1, ":", position_que + parameter[0], "\t",
                  query_seq[position_que:position_que + parameter[0]])
            position_ref += parameter[0]
            position_que += parameter[0]
        elif parameter[1] == 2:
            print("Ref:", position_ref + 1, ":", position_ref + parameter[0], "\t",
                  reference[position_ref:position_ref + parameter[0]])
            print("Que:", position_que + 1, ":", position_que + parameter[0], "\t",
                  len(query_seq[position_que:position_que + parameter[0]]) * "-")
            position_ref += parameter[0]
        elif parameter[1] == 1:
            print("Ref:", position_ref + 1, ":", position_ref + parameter[0], "\t",
                  len(reference[position_ref:position_ref + parameter[0]]) * "-")
            print("Que:", position_que + 1, ":", position_que + parameter[0], "\t",
                  query_seq[position_que:position_que + parameter[0]])
            position_que += parameter[0]
    print ()

    #Producing a continuous Ref & Que sequences (with respective in/del gaps):
    print ("Continuous aligned sequences: Reference and Query with gaps:")
    reference_with_gaps = ""
    position_ref = 0
    print ("Reference:", "\t", "( Length of seq:", len(reference), ")")
    for parameter in hit.cigar:
        if parameter[1] == 0:
            reference_with_gaps += reference[position_ref:position_ref + parameter[0]]
            position_ref += parameter[0]
        elif parameter[1] == 2:
            reference_with_gaps += reference[position_ref:position_ref + parameter[0]]
            position_ref += parameter[0]
        elif parameter[1] == 1:
            reference_with_gaps += len(reference[position_ref:position_ref + parameter[0]]) * "-"
    print (reference_with_gaps, "\t", len(reference_with_gaps), "\t", reference_with_gaps.count("-"), "insertions")
    query_with_gaps = ""
    position_que = 0
    print ("Query Seq:", "\t", "( Length of seq:", len(query_seq), ")")
    for parameter in hit.cigar:
        if parameter[1] == 0:
            query_with_gaps += query_seq[position_que:position_que + parameter[0]]
            position_que += parameter[0]
        elif parameter[1] == 2:
            query_with_gaps += len(query_seq[position_que:position_que + parameter[0]]) * "-"
        elif parameter[1] == 1:
            query_with_gaps += query_seq[position_que:position_que + parameter[0]]
            position_que += parameter[0]
    print (query_with_gaps, "\t", len(query_with_gaps), "\t", query_with_gaps.count("-"), "deletions")
    print ()

    # Summarizing the total number of matches/mismatches, insertions and deletions:
    mis_matches = 0
    insertions = 0
    deletions = 0
    for k in hit.cigar:
        if k[1] == 0:
            mis_matches += k[0]
        elif k[1] == 1:
            insertions += k[0]
        elif k[1] == 2:
            deletions += k[0]
    print ("Number of Matches/Mismatches:", mis_matches, "\t", "Number of Insertions:", insertions, "\t", "Number of deletions:", deletions)
    print ()

    # Prints gaps in the reference - i.e. insertions are included, this is not correct!
    print("R e f e r e n c e  w i t h  g a p s:")
    position_on_contig = hit.r_st
    for base_ref, base_que in zip(reference_with_gaps, query_with_gaps):
        print("Position on contig:", position_on_contig, "\t", "Ref Base:", base_ref, "\t", "Que Base:", base_que)
        position_on_contig += 1

    #Printing out each individual base in ref and seq on contig position:
        #(prints no gaps in the reference; modified for insertions)
    print("Dissection of each base with position on the contig:")
    position_on_contig = hit.r_st
    for base_ref, base_que in zip(reference_with_gaps, query_with_gaps):
        if base_ref != "-":
            print("Position on contig:", position_on_contig, "\t", "Ref Base:", base_ref, "\t", "Que Base:", base_que)
            position_on_contig += 1
    print ()

    # Match/mismatches check: in the aligned sequences:
    print("Summary of mismatches:")
    position_contig_mismatch = hit.r_st
    mismatch_counter = 0
    for base_ref, base_que in zip(reference_with_gaps, query_with_gaps):
        if base_ref != "-" and base_que != "-":
            if base_ref != base_que:
                print("Mismatch at contig position:", position_contig_mismatch, "\t", base_ref, "->", base_que)
                mismatch_counter += 1
        elif base_ref == "-":
            position_contig_mismatch -= 1
        position_contig_mismatch += 1
    print ("\t", "Total mismatches:", mismatch_counter)
    print ()
    return hit, reference, query_seq, reference_with_gaps, query_with_gaps


    # P O S I T I V E  S T R A N D S :
#AlignQueToRef("/Users/kristinaulicna/Documents/Rotation_1_Rob_Lowe/Nanopore sequences(fast5)/CLN_SMD_048987_20180430_FAH86090_MN27963_sequencing_run_BeeRapidNormal_81579_read_8_ch_239_strand.fast5")
#AlignQueToRef("/Users/kristinaulicna/Documents/Rotation_1_Rob_Lowe/Nanopore sequences(fast5)/CLN_SMD_048987_20180430_FAH86090_MN27963_sequencing_run_BeeRapidNormal_81579_read_6218_ch_429_strand.fast5")

    # N E G A T I V E  S T R A N D S :
#AlignQueToRef("/Users/kristinaulicna/Documents/Rotation_1_Rob_Lowe/Nanopore sequences(fast5)/CLN_SMD_048987_20180430_FAH86090_MN27963_sequencing_run_BeeRapidNormal_81579_read_16_ch_305_strand.fast5")
#AlignQueToRef("/Users/kristinaulicna/Documents/Rotation_1_Rob_Lowe/Nanopore sequences(fast5)/CLN_SMD_048987_20180430_FAH86090_MN27963_sequencing_run_BeeRapidNormal_81579_read_15_ch_151_strand.fast5")

help(AlignQueToRef)



# Calculation of the 'standardisation factor':

import h5py
file = h5py.File("/Users/kristinaulicna/Documents/Rotation_1_Rob_Lowe/Nanopore sequences(fast5)/CLN_SMD_048987_20180430_FAH86090_MN27963_sequencing_run_BeeRapidNormal_81579_read_6_ch_479_strand.fast5", "r")

import numpy as np
readSeq = np.array(file["Analyses"]["Basecall_1D_000"]["BaseCalled_template"]["Fastq"])
readSeq = str(readSeq)
readSeq = readSeq.split("\\n")
readSeq = readSeq[1]
print (type(readSeq), "\t", len(readSeq), "\t", readSeq)

print("ReadSeq 1-120", "\t", readSeq[0:120])

signal = file["Raw"]["Reads"]["Read_6"]["Signal"]
print (signal)
print (signal[0:600])       #because I need 5 data points for 1 base in range [0:120]

events = np.array(file["Analyses"]["Basecall_1D_000"]["BaseCalled_template"]["Events"])
print("Events in total:", len(events))
#print (events[0:1000])       #range picked randomly, you don't know how many repetitive values there are per single base...


factor_list = list()
signal_order = 0
events_order = 0
while signal_order < len(signal) and events_order < len(events):
    factor = (np.mean(signal[signal_order:signal_order+5])) / events[events_order][0]
    factor_list.append(factor)
    signal_order += 5
    events_order += 1
mean_factor_value = np.mean(factor_list[0:len(factor_list)])
stdev_factor_value = np.std(factor_list[0:len(factor_list)])
print ()
print ("Mean Factor Value:", mean_factor_value, "\t", "St.Dev.:", stdev_factor_value)