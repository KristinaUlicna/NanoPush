import h5py
import numpy as np
import mappy as mp
import os
import time


# All k-mer combinations for dictionary keys:
bases = ["A", "C", "G", "T"]
bases = [bases for i in range(5)]

combos = []
for base_5 in bases[4]:
    for base_4 in bases[3]:
        for base_3 in bases[2]:
            for base_2 in bases[1]:
                for base_1 in bases[0]:
                    combos.append(base_1 + base_2 + base_3 + base_4 + base_5)
print ("\nCombinations for keys:", len(combos), combos)


# Create the Function:

print("Loading reference genome...\n")
honeybee_genome = mp.Aligner("/Users/kristinaulicna/Documents/Rotation_1/GCF_000002195.4_Amel_4.5_genomic.fa.gz",
                           preset="map-ont")

new_Nanopolish_model_dict = {}

def NewKmerModel(fast5_file):
    """Creates a new 'Nanopolish' model as a dictionary from events data available for each k-mer."""
    input_file = h5py.File(fast5_file, "r")
    events = np.array(input_file["Analyses"]["Basecall_1D_000"]["BaseCalled_template"]["Events"])
    input_file.close()
    for line in events:
        kmer = line[4].decode('UTF-8')
        if kmer not in new_Nanopolish_model_dict:
            new_Nanopolish_model_dict[kmer] = []
        new_Nanopolish_model_dict[kmer].append(float(line[0]))
    return new_Nanopolish_model_dict


# Calling the function by iterating over 3,526 files on my Mac (Worker_1):

start_time = time.process_time()

dr = "/Users/kristinaulicna/Documents/Rotation_1/NanoporeSeq_fast5/"
folder = os.listdir(dr)

counter = 0
for fastQ_file in folder:
    dictionary = NewKmerModel(dr + fastQ_file)
    #print("File:", fastQ_file)
    counter += 1

print ("\nDictionary keys:", len(dictionary.keys()))
print ("Dictionary values:", len(dictionary.values()))

coverage_per_kmer = []
over_100_counter = 0
for key, value in dictionary.items():
    coverage_per_kmer.append(len(value))
    if len(value) >= 100:
        over_100_counter += 1

print ("\nNumber of k-mers covered in the dictionary:", len(dictionary.keys()))
print ("Total:", 4**5, ", therefore", 4**5 - len(dictionary.keys()), "kmer were not covered.")
print ("Coverage of values per kmer:", "\t", coverage_per_kmer)
print ("Number of kmers with more than 100 datapoints:", "\t", over_100_counter)

for combo in combos:
    if combo not in dictionary.keys():
        print ("These kmers not covered:", combo)
print ("All combinations for kmers checked.")
        #Proof: print (dictionary[combo])


print ("\nTotal number of FAST5 files:", counter)
print ("Processing took", time.process_time() - start_time, "seconds, i.e.", (time.process_time() - start_time) / 60, "minutes.")
print ("Processing rate is", ((time.process_time() - start_time) / counter), "seconds per file, or", (counter / (time.process_time() - start_time)), "files per second.")


# Average the values per each kmer & write a file: [kmer] : [mean_signal, stdev_signal]

for key, value in dictionary.items():
    print ("Example item of a dictionary: BEFORE", key, value)
    break

for key, value in dictionary.items():
    dictionary[key] = [np.mean(value), np.std(value)]

for key, value in dictionary.items():
    print ("Example item of a dictionary: AFTER", key, value)
    break


# Write a file to store these data:

file = open("/Users/kristinaulicna/Documents/Rotation_1/NanoPush_5-mer_Model.txt", "w")
for key, value in dictionary.items():
    string = ''
    string += (str(key) + "\t" + str(value[0]) + "\t" + str(value[1]) + "\n")
    file.write(string)
file.close()

# The output file is not ordered!
    # Sort (=alphabetically order) in terminal:
    # sort New_Model_Kmer_Values.txt -k 1,1 > Done.txt