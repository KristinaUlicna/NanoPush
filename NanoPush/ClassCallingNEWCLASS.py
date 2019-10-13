# Check if AlignQueToRef class works on all types of strands:

from NEWCLASS import Alignment
import os
import time

start_time = time.process_time()

dr = "/Users/kristinaulicna/Documents/Rotation_1/NanoporeSequences_fast5/"
folder = os.listdir(dr)

faulty_files = []
faulty_counter = 0
counter = 0
for fastQ_file in folder:
    call = Alignment(dr + fastQ_file, "/Users/kristinaulicna/Documents/Rotation_1/NanoporeSequencesCpGsNEWnanopolish_txt/")
    print("File:", fastQ_file, call.readID)
    try:
        call.AlignQueToRef()
    except:
        fastQ_file = fastQ_file.split("_")
        faulty_files.append(fastQ_file[10] + "_" + fastQ_file[11] + "_" + fastQ_file[12] + "_" + fastQ_file[13])
        faulty_counter += 1
        continue
    try:
        call.ExtractCpGsignal()
    except:
        fastQ_file = fastQ_file.split("_")
        faulty_files.append(fastQ_file[10] + "_" + fastQ_file[11] + "_" + fastQ_file[12] + "_" + fastQ_file[13])
        faulty_counter += 1
        continue
    call.StoreFileOutput()
    counter += 1
print ("\nTotal number of FAST5 files:", counter+faulty_counter)
print ("\nNumber of files analysed:", counter)
print ("\nNumber of faulty files:", faulty_counter)
print ("\nList of faulty files:", faulty_files)
print ("\nProcessing took", time.process_time() - start_time, "seconds, i.e.", (time.process_time() - start_time) / 60, "minutes.")
print ("\nProcessing rate is", ((time.process_time() - start_time) / counter+faulty_counter), "seconds per file, or", (counter+faulty_counter / (time.process_time() - start_time)), "files per second.")

# For 'cat | sort > .txt' processing in Terminal, files need to be stripped of any extra lines:
    # File name, header, descriptions, etc. - only rows with columns should stay!

"""
    # P O S I T I V E   S T R A N D:
c1 = Alignment("/Users/kristinaulicna/Documents/Rotation_1/NanoporeSequences_fast5/CLN_SMD_048987_20180430_FAH86090_MN27963_sequencing_run_BeeRapidNormal_81579_read_427_ch_225_strand.fast5", "/Users/kristinaulicna/Documents/Rotation_1/")
call_1a = c1.AlignQueToRef()
#for line in call_1a:
#    print (line)
call_1b = c1.ExtractCpGsignal()
#for line in call_1b:
#    print (line)
call_1c = c1.StoreFileOutput()


    # N E G A T I V E   S T R A N D:
c2 = Alignment("/Users/kristinaulicna/Documents/Rotation_1_Rob_Lowe/Nanopore sequences(fast5)/CLN_SMD_048987_20180430_FAH86090_MN27963_sequencing_run_BeeRapidNormal_81579_read_25_ch_177_strand.fast5")
call_2a = c2.AlignQueToRef()
call_2b = c2.ExtractCpGsignal()
call_2c = c2.StoreFileOutput()


    # M U L T I P L E   H I T S:
c3 = Alignment("/Users/kristinaulicna/Documents/Rotation_1_Rob_Lowe/Nanopore sequences(fast5)/CLN_SMD_048987_20180430_FAH86090_MN27963_sequencing_run_BeeRapidNormal_81579_read_8_ch_239_strand.fast5")
call_3a = c3.AlignQueToRef()
call_3b = c3.ExtractCpGsignal()
call_3c = c3.StoreFileOutput()

"""