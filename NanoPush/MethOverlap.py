# TODO: Create a list of 100 CpGs which are methylated and 100 CpGs which are unmethylated.

# Whole-genome bisulfide sequencing data (create a dictionary, not a list!!!):

with open("/Users/kristinaulicna/Documents/Rotation_1/GSM1360877_worker_brain_10x_meth.txt", "r") as WGBSeq_data:
    WGBSeq_list = []
    for line in WGBSeq_data:
        WGBSeq_list.append(line.strip().split('\t'))
WGBSeq_list = WGBSeq_list[::-1]

print ("WGBSeq data: [Contig, Position, Coverage, Methylation_status]:")
for counter, mini_list in enumerate(WGBSeq_list):
    print (mini_list)
print ("Counter:", counter)

"""
# Maxifile (cat&sorted) with drone data:

sorted_drone = []
for counter, CpG_line in enumerate(open("/Users/kristinaulicna/Documents/Rotation_1/Sorted_Drone_CpGs.txt", "r")):
    if counter <= 500000:
        continue
    sorted_drone.append(CpG_line.strip().split('\t'))
    if counter >= 1000000:
        break
print ("First line:", sorted_drone[0])
print ("Lines in sorted_drone:", len(sorted_drone))

for line in sorted_drone:
    if line[1] == '7903303':
        print (line)
"""
"""
print ("\nMaxifile sorted drone data:")
for counter, CpG_line in enumerate(sorted_drone):
    print (CpG_line)
    #if counter >= 2000:
    #    break
print ("Counter:", counter)
"""

# Select those CpGs from Maxifile which are methylated & which are unmethylated:

methylated_CpG = []
unmethylated_CpG = []
"""
while len(methylated_CpG) <= 100:
    for CpG_line in sorted_drone:
        for line in WGBSeq_list:
            if CpG_line[0] == line[0] and CpG_line[1] == line[1]:
                methylated_CpG.append(sorted_drone[CpG_line])
            print (methylated_CpG)
        break
"""
print ("\nMethylated CpG:", len(methylated_CpG), methylated_CpG)
print ("\nUnethylated CpG:", len(unmethylated_CpG), unmethylated_CpG)