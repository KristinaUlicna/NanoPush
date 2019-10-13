# General directory for HDF5 files:
"/Users/kristinaulicna/Documents/Rotation_1_Rob_Lowe/Nanopore sequences(fast5)/..."


# UNIVERSAL FUNCTION FOR THE ENTIRE ANALYSIS OF ALL (+ve or -ve) HDF5 FILES:

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
    if hit.strand == -1:
        query = readSeq[::-1]
        query = str(query) \
            .replace('T', '%temp1%').replace('A', 'T').replace('%temp1%', 'A') \
            .replace('C', '%temp2%').replace('G', 'C').replace('%temp2%', 'G')
        q_st = len(readSeq) - hit.q_en
        q_en = len(readSeq) - hit.q_st
        print ("Strand in negative orientation. Converted seq:", query)
        print ("New q_st & q_en:", q_st, q_en)

    else:
        query = readSeq
        q_st = hit.q_st
        q_en = hit.q_en
        print ("Strand in positive orientation. Unchanged seq:", query)
        print ("Old q_st & q_en:", q_st, q_en)

    reference = alignment.seq(hit.ctg, hit.r_st, hit.r_en)
    query_seq = query[q_st:q_en]
    print ("Reference:", "\t", len(reference), "\t", reference)
    print ("Query seq:", "\t", len(query_seq), "\t", query_seq)
    print ()

    #Breaking down aligned Ref & Que sequences accoring to the Cigar Parameters:
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
    print("Number of Matches/Mismatches:", mis_matches, "\t",
          "Number of Insertions:", insertions, "\t",
          "Number of deletions:", deletions)
    print ()

    #Position of each base on the contig + signal events data:
    print ("Events signal data + position of each base on the contig:", hit.ctg)
    #print ()
    #print ("CtgP " + "QueP " + "EveP " + "Signal " + "Shift " + "RefB " + "QueB " + "KmerB " + "Kmer ")
    events = np.array(file["Analyses"]["Basecall_1D_000"]["BaseCalled_template"]["Events"])
    reference_with_gaps = ""
    position_ref = 0
    query_with_gaps = ""
    position_que = 0
    position_events = 1
    disaster_counter = 0
    itrevents = iter(events)
    while position_events < q_st:
        line = next(itrevents)
        if line[5] != 0:
            position_events += line[5]
    for parameter in hit.cigar:
        if parameter[1] == 0:
            query_with_gaps += query_seq[position_que:position_que + parameter[0]]
            reference_with_gaps += reference[position_ref:position_ref + parameter[0]]
            for i in range(parameter[0]):
                refBase = reference[position_ref]
                queBase = query_seq[position_que]
                if refBase == queBase:
                    while position_events <= (position_que + q_st):
                        if position_events == (position_que + q_st):
                            k_mer = str(line[4])
                            #print(str(hit.r_st + position_ref) + "\t" + str(q_st + position_que) + "\t"
                            #      + str(position_events) + "\t" + str(line[0]) + "\t" + str(line[5]) + "\t"
                            #      + str(refBase) + "\t" + str(queBase) + "\t" + k_mer[4] + "\t" + k_mer[2:7])
                            #if q_st + position_que >= 10000:
                            #    break
                            if refBase != k_mer[4] and queBase != k_mer[4] :
                                disaster_counter += 1
                        line = next(itrevents)
                        if line[5] != 0:
                            position_events += line[5]
                position_ref += 1
                position_que += 1
        elif parameter[1] == 2:
            query_with_gaps += len(query_seq[position_que:position_que + parameter[0]]) * "-"
            reference_with_gaps += reference[position_ref:position_ref + parameter[0]]
            position_ref += parameter[0]
        elif parameter[1] == 1:
            query_with_gaps += query_seq[position_que:position_que + parameter[0]]
            position_que += parameter[0]
            reference_with_gaps += len(reference[position_ref:position_ref + parameter[0]]) * "-"
    print ("Disaster counter: ", disaster_counter)
    print ()

    #Producing continuous Ref & Que sequences (with respective in/del gaps):
    print ("Continuous aligned sequences: Reference and Query with gaps:")
    print ("Reference with gaps:", "\t", len(reference_with_gaps), "\t", reference_with_gaps.count("-"), "ins's", "\t", reference_with_gaps)
    print ("Query Seq with gaps:", "\t", len(query_with_gaps), "\t", query_with_gaps.count("-"), "del's", "\t", query_with_gaps)
    print ()

    #Match/mismatches check in the aligned sequences:
    print("Summary of mismatches:")
    position_contig_mismatch = hit.r_st
    mismatch_counter = 0
    for base_ref, base_que in zip(reference_with_gaps, query_with_gaps):
        if base_ref != "-" and base_que != "-":
            if base_ref != base_que:
                #print("Mismatch at contig position:", position_contig_mismatch, "\t", base_ref, "->", base_que)
                mismatch_counter += 1
        elif base_ref == "-":
            position_contig_mismatch -= 1
        position_contig_mismatch += 1
    print ("\t", "Total mismatches:", mismatch_counter)
    print ()


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

# It is more efficient to open the huge file with the entire honeybee genome once and then just align the sequences to it,
    # rather than opening the file everytime the new strand is being mapped:
    # therefore, put it ourside of your class because it is not file-specific handling & so doesn't need to be repeated


import mappy as mp
print ("Loading genome reference...\n")
alignment = mp.Aligner("/Users/kristinaulicna/Documents/Rotation_1/GCF_000002195.4_Amel_4.5_genomic.fa.gz",
                           preset="map-ont")


    # Problematic file (drone):
AlignQueToRef("/Users/kristinaulicna/Documents/Rotation_1/ProblematicRead/CLN_SMD_048987_20180828_FAJ16231_MN27963_sequencing_run_bee_drone_ampure_spike_72529_read_44800_ch_343_strand.fast5", alignment)

    # Control file (worker, positive):
AlignQueToRef("/Users/kristinaulicna/Documents/Rotation_1/NanoporeSequences_fast5/CLN_SMD_048987_20180430_FAH86090_MN27963_sequencing_run_BeeRapidNormal_81579_read_6_ch_479_strand.fast5", alignment)

    # Control file (worker, negative):
AlignQueToRef("/Users/kristinaulicna/Documents/Rotation_1/NanoporeSequences_fast5/CLN_SMD_048987_20180430_FAH86090_MN27963_sequencing_run_BeeRapidNormal_81579_read_24_ch_347_strand.fast5", alignment)

