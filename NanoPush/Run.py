    #Choose which Alignment you want (old = AlignQueToRef, new = AlignClassNew):
from AlignQueToRef import Alignment
from AlignClassNew import Alignment
import sys
import os
import time

input_folder = sys.argv[1]
output_folder = sys.argv[2]

folder = os.listdir(input_folder)
start_time = time.process_time()
faulty_files = []
faulty_counter = 0
counter = 0
for fastQ_file in folder:
    call = Alignment(input_folder + "/" + fastQ_file, output_folder)
    print(fastQ_file, call.readID)
    try:
        call.AlignQueToRef()
    except:
        fastQ_file = fastQ_file.split("_")
        #faulty_files.append(fastQ_file[10] + "_" + fastQ_file[11] + "_" + fastQ_file[12] + "_" + fastQ_file[13])
        faulty_counter += 1
        continue
    try:
        call.ExtractCpGsignal()
    except:
        fastQ_file = fastQ_file.split("_")
        #faulty_files.append(fastQ_file[10] + "_" + fastQ_file[11] + "_" + fastQ_file[12] + "_" + fastQ_file[13])
        faulty_counter += 1
        continue
    call.StoreFileOutput()
    counter += 1
print ("\nTotal number of FAST5 files:", counter+faulty_counter)
print ("\nNumber of files analysed:", counter)
print ("\nNumber of faulty files:", faulty_counter)
print ("\nProcessing took", time.clock() - start_time, "seconds, i.e.", (time.clock() - start_time) / 60, "minutes.")
print ("\nProcessing rate is", ((time.clock() - start_time) / counter+faulty_counter), "seconds per file, or", (counter+faulty_counter / (time.clock() - start_time)), "files per second.")



