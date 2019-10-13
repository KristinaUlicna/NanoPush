from AlignQueToRef import Alignment

call = Alignment("/Users/kristinaulicna/Documents/Rotation_1/CLN_SMD_048987_20180828_FAJ16231_MN27963_sequencing_run_bee_drone_ampure_spike_72529_read_142623_ch_205_strand.fast5", "/Users/kristinaulicna/Documents/Rotation_1/")

for line in call.AlignQueToRef():
    if line[0][0:5] == "56504":
        print (line)

for line in call.ExtractCpGsignal():
    print (line)