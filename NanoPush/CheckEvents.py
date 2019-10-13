# ----------- MODIFY THE FUNCTION FOR THE NEGATIVE STRAND:

def AlignQueToRef(HDF5_file, alignment):
    """Insert HDF5 file name into "" with a full directory (from /Users/...)
    1.) Aligns the FastQ sequence (ReadSeq) to the reference sequence (honeybee genome).
    2.) Reverse-complements strand if in negative orientation & changes q_st and q_en.
    3.) Prints the reference and query sequences.
    4.) Chops the sequence according to the cigar parameters (mism, in/dels + summary).
    5.) Produces a table with contig positions and base calls for each match event recorded:
        (format: CtgPos, QuePos, EvePos, Signal, Shift, RefBase, QueBase, KmerBase, Kmer)
    6.) Prints the reference and query sequences with gaps (modified for in/dels).
    7.) Produces a summary of mismatches & their total count."""
    print ("HDF5 file name:", "\t", HDF5_file)

    #Extraction of FastQ sequence (string) from an HDF5 file type:
    import h5py
    file = h5py.File(HDF5_file, "r")
    import numpy as np
    readSeq = np.array(file["Analyses"]["Basecall_1D_000"]["BaseCalled_template"]["Fastq"])
    readSeq = str(readSeq)
    readSeq = readSeq.split("\\n")
    readSeq = readSeq[1]
    print ("ReadSeq:", "\t", len(readSeq), "\t", readSeq)

    #Alignment of the readSeq sequence to the reference genome (honeybee - Apis mellifera):
    for hit in alignment.map(readSeq):
        print ("Hit parameters:", "\t", hit)

    #Strand conversion into positive orientation + hit.q_st & hit.q_en update into q_st & q_en:
    reference = alignment.seq(hit.ctg, hit.r_st, hit.r_en)
    if hit.strand == -1:
        reference = reference[::-1]
        reference = str(reference) \
            .replace('T', '%temp1%').replace('A', 'T').replace('%temp1%', 'A') \
            .replace('C', '%temp2%').replace('G', 'C').replace('%temp2%', 'G')
    query_seq = readSeq[hit.q_st:hit.q_en]
    print ("Reference:", "\t", len(reference), "\t", reference)
    print ("Query seq:", "\t", len(query_seq), "\t", query_seq)
    print ()

    #Position of each base on the contig + signal events data:
    print ("Events signal data + position of each base on the contig:", hit.ctg)
    print ()
    print ("CtgP " + "QueP " + "EveP " + "Signal " + "Shift " + "RefB " + "QueB " + "KmerB " + "Kmer " + ":")
    events = np.array(file["Analyses"]["Basecall_1D_000"]["BaseCalled_template"]["Events"])
    my_list = []
    position_ref = 0
    position_que = 0
    position_events = 1
    itrevents = iter(events)
    while position_events < hit.q_st:
        line = next(itrevents)
        if line[5] != 0:
            position_events += line[5]
    for parameter in hit.cigar:
        if parameter[1] == 0:
            for i in range(parameter[0]):
                refBase = reference[position_ref]
                queBase = query_seq[position_que]
                if refBase == queBase:
                    while position_events <= (position_que + hit.q_st):
                        if position_events == (position_que + hit.q_st):
                            k_mer = str(line[4])
                            my_list.append([str(hit.r_st + position_ref), str(hit.q_st + position_que),
                                            str(position_events), str(line[0]), str(line[5]),
                                            str(refBase), str(queBase), k_mer[4], k_mer[2:7]])
                            #my_list.append("\n")
                            if refBase != k_mer[4] and queBase != k_mer[4]:
                                raise Exception("MisMatch", "RefBase, QueBase and K-merBase are not matching",
                                                position_events)
                        line = next(itrevents)
                        if line[5] != 0:
                            position_events += line[5]
                position_ref += 1
                position_que += 1
        elif parameter[1] == 2:
            position_ref += parameter[0]
        elif parameter[1] == 1:
            position_que += parameter[0]
    print (my_list)
    return my_list



#TODO: How do I pull out data for an individual contig position? Say, I am interested in position xx,xxx,xxx - find out!
    # Not possible to do from printed out data. You could store the output in a dictionary, such as in Rob's example:
    # (written in terminal, so different structure)
"""    
>>>dict = {3636:["A","C"], 3637:["C","C"]}
>>>dict[3636]
['A', 'C']
>>> dict[3637]
['C', 'C']
>>> dict={}
>>> dict[3636]=["A","C"]
>>> dict[3637]
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
KeyError: 3637
>>> dict[3636]
['A', 'C']
>>> lst=[[3636,"A","C"],[3637,"C","C"]]
>>> lst[0][0]
3636
>>> 3637-lst[0][0]
1
>>> lst[3637-lst[0][0]]
[3637, 'C', 'C']
>>> lst=[[3636,"A","C"],[3637,"C","C"],[3638,"A","A"]]
>>> lst[3637-lst[0][0]]
[3637, 'C', 'C']
>>> lst[3638-lst[0][0]]
[3638, 'A', 'A']
"""

    # P O S I T I V E  S T R A N D S :
#AlignQueToRef("/Users/kristinaulicna/Documents/Rotation_1_Rob_Lowe/Nanopore sequences(fast5)/CLN_SMD_048987_20180430_FAH86090_MN27963_sequencing_run_BeeRapidNormal_81579_read_6_ch_479_strand.fast5")
#AlignQueToRef("/Users/kristinaulicna/Documents/Rotation_1_Rob_Lowe/Nanopore sequences(fast5)/CLN_SMD_048987_20180430_FAH86090_MN27963_sequencing_run_BeeRapidNormal_81579_read_6218_ch_429_strand.fast5")

    # N E G A T I V E  S T R A N D S :
#AlignQueToRef("/Users/kristinaulicna/Documents/Rotation_1_Rob_Lowe/Nanopore sequences(fast5)/CLN_SMD_048987_20180430_FAH86090_MN27963_sequencing_run_BeeRapidNormal_81579_read_16_ch_305_strand.fast5")
#AlignQueToRef("/Users/kristinaulicna/Documents/Rotation_1_Rob_Lowe/Nanopore sequences(fast5)/CLN_SMD_048987_20180430_FAH86090_MN27963_sequencing_run_BeeRapidNormal_81579_read_15_ch_151_strand.fast5")


# It is more efficient to open the huge file with the entire honeybee genome once and then just align the sequences to it,
    # rather than opening the file everytime the new strand is being mapped:
    # therefore, put it ourside of your class because it is not file-specific handling & so doesn't need to be repeated


import mappy as mp
print ("Loading genome reference...")
alignment = mp.Aligner("/Users/kristinaulicna/Documents/Rotation_1_Rob_Lowe/GCF_000002195.4_Amel_4.5_genomic.fa.gz",
                           preset="map-ont")
print ("Reference genome loaded.")
print ()


    # Positive strands:
AlignQueToRef("/Users/kristinaulicna/Documents/Rotation_1_Rob_Lowe/Nanopore sequences(fast5)/CLN_SMD_048987_20180430_FAH86090_MN27963_sequencing_run_BeeRapidNormal_81579_read_28_ch_21_strand.fast5", alignment)

    # Negative strands:
#AlignQueToRef("/Users/kristinaulicna/Documents/Rotation_1_Rob_Lowe/Nanopore sequences(fast5)/CLN_SMD_048987_20180430_FAH86090_MN27963_sequencing_run_BeeRapidNormal_81579_read_24_ch_461_strand.fast5", alignment)
#AlignQueToRef("/Users/kristinaulicna/Documents/Rotation_1_Rob_Lowe/Nanopore sequences(fast5)/CLN_SMD_048987_20180430_FAH86090_MN27963_sequencing_run_BeeRapidNormal_81579_read_25_ch_177_strand.fast5", alignment)
AlignQueToRef("/Users/kristinaulicna/Documents/Rotation_1_Rob_Lowe/Nanopore sequences(fast5)/CLN_SMD_048987_20180430_FAH86090_MN27963_sequencing_run_BeeRapidNormal_81579_read_26_ch_215_strand.fast5", alignment)


"""
AlignQueToRef("/Users/kristinaulicna/Documents/Rotation_1_Rob_Lowe/Nanopore sequences(fast5)/CLN_SMD_048987_20180430_FAH86090_MN27963_sequencing_run_BeeRapidNormal_81579_read_28_ch_235_strand.fast5", alignment)
AlignQueToRef("/Users/kristinaulicna/Documents/Rotation_1_Rob_Lowe/Nanopore sequences(fast5)/CLN_SMD_048987_20180430_FAH86090_MN27963_sequencing_run_BeeRapidNormal_81579_read_28_ch_269_strand.fast5", alignment)
AlignQueToRef("/Users/kristinaulicna/Documents/Rotation_1_Rob_Lowe/Nanopore sequences(fast5)/CLN_SMD_048987_20180430_FAH86090_MN27963_sequencing_run_BeeRapidNormal_81579_read_29_ch_313_strand.fast5", alignment)
AlignQueToRef("/Users/kristinaulicna/Documents/Rotation_1_Rob_Lowe/Nanopore sequences(fast5)/CLN_SMD_048987_20180430_FAH86090_MN27963_sequencing_run_BeeRapidNormal_81579_read_30_ch_271_strand.fast5", alignment)
AlignQueToRef("/Users/kristinaulicna/Documents/Rotation_1_Rob_Lowe/Nanopore sequences(fast5)/CLN_SMD_048987_20180430_FAH86090_MN27963_sequencing_run_BeeRapidNormal_81579_read_31_ch_102_strand.fast5", alignment)
AlignQueToRef("/Users/kristinaulicna/Documents/Rotation_1_Rob_Lowe/Nanopore sequences(fast5)/CLN_SMD_048987_20180430_FAH86090_MN27963_sequencing_run_BeeRapidNormal_81579_read_31_ch_444_strand.fast5", alignment)

"""