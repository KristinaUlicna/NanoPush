#Create data set from SELECTED CpGs from Rob (200 meth & unmeth) acc. to the WGBSeq_dict:

from FileHandling_Functions import FilterNone
#FilterNone("/Users/kristinaulicna/Documents/Rotation_1/Archive/DRONE/NewSelectedFromRob/Sorted_Drone_Sel_All.txt", )
#FilterNone("/Users/kristinaulicna/Documents/Rotation_1/Archive/DRONE/NewSelectedFromRob/Sorted_Drone_Sel_M_All.txt", "/Users/kristinaulicna/Documents/Rotation_1/Archive/DRONE/NewSelectedFromRob/HelperFolder/Sorted_Drone_Sel_M_NoNone.txt")


#TODO: Before extraction, move the '_NoNone.txt' file into a separate folder (e.g. /HelperFolder/):
import time
from ExtractUnMethCpGs_Function import ExtractCpGsNewModel

start_time = time.process_time()
unmeth_list, meth_list = ExtractCpGsNewModel("/Users/kristinaulicna/Documents/Rotation_1/Archive/DRONE/NewSelectedFromRob/HelperFolder/")
print ("\nUnmeth_list:", len(unmeth_list), "\n\tMeth_List:", len(meth_list))
print ("\nProcessing took", time.process_time() - start_time, "seconds, i.e.", (time.process_time() - start_time) / 60, "minutes.")


#Create files from which ML will source:

unmeth_file = open("/Users/kristinaulicna/Documents/Rotation_1/Archive/DRONE/NewSelectedFromRob/CpGs_sel_Drone_meth_less_0_1.txt", "w")
string = ''
for counter, single_line in enumerate(unmeth_list):
    if counter < 500:
        for parameter in single_line:
            string += (str(parameter) + "\t")
        string = string[:-1]
        string += "\n"
unmeth_file.write(string)
unmeth_file.close()

meth_file = open("/Users/kristinaulicna/Documents/Rotation_1/Archive/DRONE/NewSelectedFromRob/CpGs_sel_Drone_meth_over_0_9.txt", "w")
string = ''
for counter, single_line in enumerate(meth_list):
    for parameter in single_line:
        string += (str(parameter) + "\t")
    string = string[:-1]
    string += "\n"
meth_file.write(string)
meth_file.close()
