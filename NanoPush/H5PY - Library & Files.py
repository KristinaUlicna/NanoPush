# ---------- Exploring the file structure:

# Directory to open files:
"""/Users/kristinaulicna/Documents/Rotation_1_Rob_Lowe/Nanopore sequences(fast5)/CLN... .fast5"""


import h5py
f1 = h5py.File("/Users/kristinaulicna/Documents/Rotation_1_Rob_Lowe/Nanopore sequences(fast5)/CLN_SMD_048987_20180430_FAH86090_MN27963_sequencing_run_BeeRapidNormal_81579_read_8_ch_239_strand.fast5", "r")

# h5py.File acts like a Python dictionary, thus we can check the keys:
print (list(f1.keys()))      # to list the individual keys
print (list(f1.values()))

for key in f1.keys():       # ...to iterate over them
    print(key)              # names of the groups in HDF5 file
print ()



# ---------- Indexing and slicing the dictionary - exploring individual layers:

# I can see there is a key called 'UniqueGlobalKey' - what is in it? Let's have a look:
# nomenclature of variable = key_layer_member

key_1_1 = f1['UniqueGlobalKey']
print ("\n", "key_1_1", key_1_1)             # summarizes the key characteristics, e.g number of members (here, 3)...
print (list(key_1_1))                   # prints the list of how the 3 members are called

key_1_2 = f1['PreviousReadInfo']
print ("\n", "key_1_2", key_1_2)             # has 0 members, therefore:
print (list(key_1_2))                   # this list is empty... (0 member names to print)

key_1_3 = f1['Raw']
print ("\n", "key_1_3", key_1_3)             # has 1 member
print (list(key_1_3))                   # 'Reads'

key_1_4 = f1['Analyses']
print ("\n", "key_1_4", key_1_4)             # has 4 members
print (list(key_1_4))                   # output: ['Basecall_1D_000', 'Calibration_Strand_Detection_000', 'RawGenomeCorrected_000', 'Segmentation_000']


#---------- What is stored under UniqueGlobalKey?

key_1_1_1 = list(f1["UniqueGlobalKey"]["tracking_id"].keys())
value_1_1_1 = list(f1["UniqueGlobalKey"]["tracking_id"].values())
print (key_1_1_1, value_1_1_1)      # both empty

key_1_1_2 = list(f1["UniqueGlobalKey"]["context_tags"].keys())
value_1_1_2 = list(f1["UniqueGlobalKey"]["context_tags"].values())
print (key_1_1_2, value_1_1_2)      # both empty

key_1_1_3 = list(f1["UniqueGlobalKey"]["channel_id"].keys())
value_1_1_3 = list(f1["UniqueGlobalKey"]["channel_id"].values())
print (key_1_1_3, value_1_1_3)      # both empty


#---------- What is stored under PreviousReadInfo?
    # Nothing, has 0 members...


#---------- What is stored under Raw? (has only 1 member; so expect your data of interest)

    # Two options how to do this:
        # 1.) Call the group by indexing - you know there is only 1 member, so index it as [0]

key_1_3_indexing = list(f1["Raw"].keys())[0]

        # 2.) Look at the names of the members which the "Raw" category contains, and print them as a list

key_1_3_naming = list(f1["Raw"]["Reads"])

print ("\n", key_1_3_indexing, key_1_3_naming, "\n")
        #output: Reads [Read_8] = this corresponds to the file name!

#---------- What is stored under Raw & Reads?

key_1_3_1 = list(f1["Raw"]["Reads"].keys())[0]
value_1_3_1 = list(f1["Raw"]["Reads"].values())[0]
print (key_1_3_1, "\t", value_1_3_1)

#---------- What is stored under Raw & Reads & Read_8?

key_1_3_1_1 = list(f1["Raw"]["Reads"]["Read_8"].keys())[0]
value_1_3_1_1 = list(f1["Raw"]["Reads"]["Read_8"].values())[0]
print (key_1_3_1_1, "\t", value_1_3_1_1)

#key_1_3_1_1_1 = list(f1["Raw"]["Reads"]["Read_8"]["Signal"].keys())[0]
#value_1_3_1_1_1 = list(f1["Raw"]["Reads"]["Read_8"]["Signal"].values())[0]
#   "Signal" category has no more keys in it, neither values (no further layer of dictionary!)

key_1_3_1_1_1 = f1["Raw"]["Reads"]["Read_8"]["Signal"]
    # output: <HDF5 dataset "Signal": shape (4224,), type "<i2">
        # i.e. data are stored in a vector with a 4224 rows and 0 columns (=their shape) in an 'i2' format (=their type)

my_data = key_1_3_1_1_1[0:4224]
print (my_data)
    # output: individual numbers!


#---------- What is stored under Analyses?

key_1_4 = list(f1["Analyses"].keys())
key_1_4_1 = list(f1["Analyses"]["Basecall_1D_000"].keys())
key_1_4_2 = list(f1["Analyses"]["Calibration_Strand_Detection_000"].keys())
key_1_4_3 = list(f1["Analyses"]["RawGenomeCorrected_000"].keys())
key_1_4_4 = list(f1["Analyses"]["Segmentation_000"].keys())

print ("\n", key_1_4, "\n", key_1_4_1, "\n", key_1_4_2, "\n", key_1_4_3, "\n", key_1_4_4)

key_1_4_1_1 = list(f1["Analyses"]["Basecall_1D_000"]["BaseCalled_template"].keys())
value_1_4_1_1 = list(f1["Analyses"]["Basecall_1D_000"]["BaseCalled_template"].values())

print (key_1_4_1_1, "\n", value_1_4_1_1)
        # Events: shape (820,)
        # Fastq: shape ()

#key_1_4_1_1_1 = list(f1["Analyses"]["Basecall_1D_000"]["BaseCalled_template"]["Events"].keys())
#key_1_4_1_1_2 = list(f1["Analyses"]["Basecall_1D_000"]["BaseCalled_template"]["Fastq"].keys())
#value_1_4_1_1_1 = list(f1["Analyses"]["Basecall_1D_000"]["BaseCalled_template"]["Events"].values())
#value_1_4_1_1_2 = list(f1["Analyses"]["Basecall_1D_000"]["BaseCalled_template"]["Fastq"].values())
# "Events" and "Fastq" category has no more keys in it, neither values (no further layer of dictionary!)

key_1_4_1_1_1 = f1["Analyses"]["Basecall_1D_000"]["BaseCalled_template"]["Events"]
key_1_4_1_1_2 = f1["Analyses"]["Basecall_1D_000"]["BaseCalled_template"]["Fastq"]
print (key_1_4_1_1_1, "\n", key_1_4_1_1_2)
        # stick to "Events" only because "Fastq" has no rows nor columns in 'shape'

my_events = key_1_4_1_1_1[0:3] #up to 820
print (my_events)



# VISUALISING THE DATA:

import numpy as np
import matplotlib.pyplot as plotting
plotting.plot(my_data[1000:4224])
plotting.show()


# Access the relevant data from another file in one go:

f3 = h5py.File("/Users/kristinaulicna/Documents/Rotation_1_Rob_Lowe/Nanopore sequences(fast5)/CLN_SMD_048987_20180430_FAH86090_MN27963_sequencing_run_BeeRapidNormal_81579_read_10_ch_358_strand.fast5", "r")
data_f3 = f3["Raw"]["Reads"]["Read_10"]["Signal"]
print (data_f3)
plotting.plot(data_f3[0:4184])  #all data
plotting.show()


# Extract the sequence data from the read ('Read_10')

import numpy as np
seq = np.array(f3["Analyses"]["Basecall_1D_000"]["BaseCalled_template"]["Fastq"])
seqString = str(seq)
seqString = seqString.split("\\n")
print (seqString)      # has '@' at the beginning

    # to print each line on a separate line:
        #for line in potato:
        #    print(line)

    # to print the raw sequence only:
        #print (seqString[1])

    # to print the sequence so that 1 line only displays 60 bases:

sequence = seqString[1]     # raw sequence string only
print (sequence, "\n",
       len(sequence), "\n",
       type(sequence), "\n")

position = 0                                # SCOPING of a variable:
while position+60 < len(sequence):
    print(sequence[position:position+60])
    position += 60

print(sequence[position:len(sequence)])
print (position)                            # while loop alters the value of 'position' = SCOPING
                                            # = after the while loop, it becomes 360
                                            # and not 0 as it was initially defined!
"""
What happened in my while loop? Martin's line-by-line explanation:

[ ..... 380 ..... ]         = string of 380 characters)
offset : 60                 = how many characters we want to print in each line

Iterations:
0 (first):    position 0, print, position 60
1 (second):   position 60, print, position 120
2 (third):    position 120, print, position 180
3 (fourth):   position 180, print, position 240
4 (fifth):    position 240, print, position 300
5 (sixth):    position 300, print, position 360
CANNOT RUN:   position 360 + 60 = 420 which is more than the len(sequence),
                therefore the condition is no longer true and while loop will terminate
                --- 
But what about the remainder 20 characters which will now not be printed?
PRINT THEM SEPARATELY! My variable 'position' has now been scoped! 
...If you print the value of position after the while loop has been executed, 
        it will no longer be '0' as defined before but '360'
Now, print the rest of the 'sequence' [360:380], i.e [position:len(sequence] with a separate print statement

And you are done! :-)

"""

# Getting the read ID

print (seqString[0])       # this is where read ID is stored
full_line = seqString[0]
print (type(full_line))
full_id_1 = full_line.split("@")[1]
full_id_2 = full_id_1.split("_")[0]
print ("Read ID:", full_id_2)               # simple way
print ("Read ID: {}".format(full_id_2))     # fancy way

"""
# BONUS: Martin's explanation of what a method vs function is:
# methods are functions callable on an object
# functions are not methods, but methods are functions
# example method: car.drive() - it is a verb of what the noun does

"""