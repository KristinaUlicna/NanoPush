# ---------- How to import Python external libraries, like H5PY, into PyCharm?

# 1.) In the Project Interpreter page of the project settings, select the desired Python interpreter or virtual environment.
# 2.) Click +
# 3.) In the Available Packages dialog box that opens, select the desired package from the list.
# 4.) If required, select the following check boxes:
# 5.) Specify Version: If this check box is selected, you can select the desired version from the drop-down list of available versions. By default, the latest version is taken.
# 6.) Options: If this check box is selected, you can type the options in the text field.
# 7.) Click Install Package.

# library = h5py
# file type = HDF5
# file extension = .fast5 (Oxford Nanopore-specific)



#---------- Accessing & plotting the signal data (picoamp values):
    # NB: an average of 5 data points (currents) is averaged and normalised (factor ~ 5.3) in the events.

import h5py
file = h5py.File("/Users/kristinaulicna/Documents/Rotation_1_Rob_Lowe/Nanopore sequences(fast5)/CLN_SMD_048987_20180430_FAH86090_MN27963_sequencing_run_BeeRapidNormal_81579_read_10_ch_358_strand.fast5", "r")

#Each file has the "Read_x" stored in its name - to access it, you would have to check each file and copy it from there.
    #An easier solution is to extract this "Read_x" from the name for each file to be analysed:

read_number = str(file)
read_number = read_number.split("_")
read_number = read_number[10] + "_" + read_number[11]
read_number = read_number.capitalize()
print ("Read number from the HDF5 file name: " + read_number)
print ()

access_signal = list(file["Raw"]["Reads"][read_number]["Signal"])[0:4000]
values_list = []
for value in access_signal:
    values_list.append(value)

#print ("Values list: ", "\t", values_list)
    #very long line, split with while loop:

print ("Signal data (40 datapoints / line):")
signal_value = 0
while signal_value + 40 < len(values_list):
    print (values_list[signal_value:signal_value+40])
    signal_value += 40
print (values_list[signal_value:len(values_list)])
print ()


import matplotlib.pyplot as plot
plot.plot(access_signal[2000:2500])
plot.show()



#---------- Accessing the query sequence from the FastQ format:

access_events = list(file["Analyses"]["Basecall_1D_000"]["BaseCalled_template"]["Events"])[0]
print ("Signal events:", "\t", access_events)

access_fastq = file["Analyses"]["Basecall_1D_000"]["BaseCalled_template"]["Fastq"]
print ("FastQ format:", "\t", access_fastq)
print ()

import numpy as np
sequence = np.array(file["Analyses"]["Basecall_1D_000"]["BaseCalled_template"]["Fastq"])
access_sequence = str(sequence)
access_sequence = access_sequence.split("\\n")
access_sequence = access_sequence[1]

#print ("Seq length:", len(access_sequence), "\t", access_sequence)
    #very long line, split with while loop:

print ("Sequence (60 characters / line):")
position = 0
while position + 60 < len(access_sequence):
    print (access_sequence[position:position+60])
    position += 60
print (access_sequence[position:len(access_sequence)])
print ()

print ("Read ID (unique for each read):")
access_readID = str(sequence)
access_readID = access_readID.split("@")
access_readID = access_readID[1]
access_readID = access_readID.split("_")
access_readID = access_readID[0]
print (access_readID)
print ()



#---------- Alignment of the query sequence to the reference sequence (Apis mellifera):

# Python tool: Minimap2 for alignment
# mappy = python module (can be imported)

import mappy as mp
alignment = mp.Aligner("/Users/kristinaulicna/Documents/Rotation_1_Rob_Lowe/GCF_000002195.4_Amel_4.5_genomic.fa.gz",
                       preset = "map-ont")
# .gz ending means that the file is zipped - you don't need to unzip it to open it in mappy!
print(alignment)

# METHOD = .seq() - takes 1 string argument (compulsory) & 2 numerical arguments (optional):
    # "NC_007080.3" = the 'RefSeq sequence' name of the honeybee linkage group (=chromosome)
    # 5 & 10 = start and end of the sequence read which you want to map to, or print
                # (it's like slicing but without the square brackets [] )
                # arguments are separated by commas!

# METHOD = .map() - takes 1 argument:
    # it's a generator!!!
    # "test" is a sliced out sequence (5:1000) from the reference genome on the contig we want to map to:
                # if we map a slice of a reference to a reference, it will give us a perfect match!
                # .map is a generator - will print the hit, containing the hit.X parameters


test = alignment.seq("NC_007080.3", 5, 1000)   # if the 'sliced' sequence is too short, it won't map!
print ("Test Seq:", type(test), test)
        # this is a test because we are pulling a 995 nucleotide-long sequence from a reference sequence
        # and mapping it back to the same reference sequence - this will give us perfect match! (as expected)

for hit in alignment.map(test):
    print (hit)

    #output: 0	995	+	NC_007080.3	14726556	5	1000	995	995	60	tp:A:P	ts:A:.	cg:Z:995M

    #        start of query alignment
    #           end of query alignment
    #               strand orientation (positive or negative)
    #                   contig (=chromosome) accession number (honeybee, Apis mellifera)
    #                               contig (=chromosome) total length (in basepairs)
    #                                           start position of mapping on reference sequence
    #                                               end position of mapping on reference sequence
    #                                                               quality score (60 = highest)
    #                                                                                   "cigar sequence" (see below)

query = alignment.map(access_sequence) #our query sequence which we pulled out from the "CLN..." file
print ("Length of Sequence to be mapped:", len(access_sequence))

for hit in alignment.map(access_sequence):
    print (hit)

    # output:  100	375	+	NC_007080.3	14726556	5067244	5067526	271	284	60	tp:A:P	ts:A:.	cg:Z:67M6D67M1D24M1I12M1D54M1D21M1I28M
                          # coincidentally, this is the same chromosome as before...

    # cigar sequence is different now: 67M6D67M1D24M1I12M1D54M1D21M1I28M
        # meaning there is a 67bp-long match or mismatch; 6bp-long deletion; 67bp-long match or mismatch; 1bp-long insertion; etc.

# Do the lengths of the query and reference sequences mapped to each other correspond?
    # you don't have to copy the values from the output to access them, just call them by their name:

print (hit.r_en - hit.r_st)     # r = reference, st = start, en = end
print (hit.q_en - hit.q_st)     # q = query

    # reference = 282; to calculate this manually, you need to add together the number of (MIS)MATCHES + DELETIONS
    # query = 275; to calculate this manually, you need to add together the number of (MIS)MATCHES + INSERTIONS


# My cigar sequence told me that the first 67 bases have been aligned one to one:
    # Are those matches or mismatches?

queSeq = access_sequence[hit.q_st:hit.q_st+67]
refSeq = alignment.seq("NC_007080.3", hit.r_st, hit.r_st+67)
print (refSeq is queSeq)
    # output: False => there is at least 1 mismatch! (position 13, G -> A)
    # check where it is:

position = 0
mismatch_counter = 0
for baseRef, baseQue in zip(refSeq, queSeq):
    position += 1
    if baseRef != baseQue:
        print ("Position:", position, "\t", "RefBase:", baseRef, "\t", "QueBase:", baseQue)
        mismatch_counter += 1
print ("Total mismatches:", mismatch_counter)

print ("Contig:", hit.ctg)     #output: 'NC_007080.3'  (RefSeq sequence name of the chromosome)
print (hit.cigar)              #output is a list of the lists,
    # where the 1st number of each sublist is the quantity (=length in bp) of the parameter; and
    # the 2nd number is the code for matches/mismatches (=0), insertions (=1) or deletions (=2)
    # why? because numbers (0, 1, 2) are computationally less demanding to store compared to characters ("M", "I", "D")



#---------- Summary of the number of matches/mismatches & in/dels in your aligned sequences:

mis_matches = 0
insertions = 0
deletions = 0
                                # these values will be scoped by the for-loop, therefore you can print them outside of the loop with their altered values
for k in hit.cigar:             # k = each individual list of the long list (remember, hit.cigar = list of lists)
    #print(k)                   # re-prints the list of the lists which we saw earlier in 'hit.cigar'
    if k[1] == 0:                   # meaning: "if the second value of each list is equal to 0..."
        mis_matches += k[0]         # meaning: "then add the first value (i.e. length of the stretch) to the overall number of mismatches
    elif k[1] == 1:
        insertions += k[0]
    elif k[1] == 2:
        deletions += k[0]

print ("Number of Matches/Mismatches:", mis_matches, "\t",
       "Number of Insertions:", insertions, "\t",
       "Number of deletions:", deletions)

print ("Ref Sequence length:", "\t", "Matches/mismatches + deletions, i.e", mis_matches, "+", deletions, "=", mis_matches+deletions, "bp", "\t", (mis_matches + deletions) == (hit.r_en - hit.r_st))
print ("Que Sequence length:", "\t", "Matches/mismatches + insertions, i.e", mis_matches, "+", insertions, "=", mis_matches+insertions, "bp", "\t", (mis_matches + insertions) == (hit.q_en - hit.q_st))

# Now, make sure your calculator works:
    # Ref seq length      = Matches/mismatches + deletions    = 273 + 9 = 282bp
    # Que seq length      = Matches/mismatches + insertions   = 273 + 2 = 275bp
            # Perfect!
