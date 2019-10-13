import time
from WGBSeqClassif_Functions import ClassifyCpGsWGBSeq

start_time = time.process_time()
ClassifyCpGsWGBSeq("/Users/kristinaulicna/Documents/Rotation_1/Archive/WORKER/Sorted_Worker_NoNone_Shuffled.txt", "/Users/kristinaulicna/Documents/Rotation_1/Archive/WORKER/", 500)
print ("File sorting took", time.process_time() - start_time, "seconds, i.e.", (time.process_time() - start_time) / 60, "minutes.")
