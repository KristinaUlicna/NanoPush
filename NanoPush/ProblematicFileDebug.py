from AlignQueToRef import Alignment

call_1 = Alignment("/Users/kristinaulicna/Documents/Rotation_1/ProblematicRead/CLN_SMD_048987_20180828_FAJ16231_MN27963_sequencing_run_bee_drone_ampure_spike_72529_read_44800_ch_343_strand.fast5", "/Users/kristinaulicna/Documents/Rotation_1/ProblematicRead/")

print ("\nCall_1.AlignQueToRef()")
for line in call_1.AlignQueToRef():
    print (line)
    if int(line[0]) >= 300:
        break

print ("\nCall_1.ExtractCpGsignal()")
for line in call_1.ExtractCpGsignal():
    print (line)



call_2 = Alignment("/Users/kristinaulicna/Documents/Rotation_1/NanoporeSequences_fast5/CLN_SMD_048987_20180430_FAH86090_MN27963_sequencing_run_BeeRapidNormal_81579_read_6_ch_479_strand.fast5", "/Users/kristinaulicna/Documents/Rotation_1/ProblematicRead/")

print ("\nCall_2.AlignQueToRef()")
for line in call_2.AlignQueToRef():
    print (line)
    if int(line[0]) >= 8394060:
        break

print ("\nCall_2.ExtractCpGsignal()")
for line in call_2.ExtractCpGsignal():
    print (line)



call_3 = Alignment("/Users/kristinaulicna/Documents/Rotation_1/NanoporeSequences_fast5/CLN_SMD_048987_20180430_FAH86090_MN27963_sequencing_run_BeeRapidNormal_81579_read_24_ch_347_strand.fast5", "/Users/kristinaulicna/Documents/Rotation_1/ProblematicRead/")

print ("\nCall_3.AlignQueToRef()")
for line in call_3.AlignQueToRef():
    print (line)
    if int(line[0]) >= 1513050:
        break

print ("\nCall_3.ExtractCpGsignal()")
for line in call_3.ExtractCpGsignal():
    print (line)