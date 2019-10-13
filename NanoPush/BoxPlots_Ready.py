#---------- Produce data sets by importing all CpGs from designated files into 2 lists:

from FileHandling_Functions import CreateDataSet
unmeth_list_worker_1 = CreateDataSet("/Users/kristinaulicna/Documents/Rotation_1/Data/CpG_all_Worker_1_meth_less_0_1.txt", 1, 500)
meth_list_worker_1 = CreateDataSet("/Users/kristinaulicna/Documents/Rotation_1/Data/CpG_all_Worker_1_meth_over_0_9.txt", 1, 272)

print("\nWORKER_1 Unmethylated List:", "\t", len(unmeth_list_worker_1), "\t", unmeth_list_worker_1[0])
print("\nWORKER_1 Methylated List:", "\t", len(meth_list_worker_1), "\t", meth_list_worker_1[0])

unmeth_list_worker_2 = CreateDataSet("/Users/kristinaulicna/Documents/Rotation_1/Data/CpG_all_Worker_2_meth_less_0_1.txt", 1, 500)
meth_list_worker_2 = CreateDataSet("/Users/kristinaulicna/Documents/Rotation_1/Data/CpG_all_Worker_2_meth_over_0_9.txt", 1, 354)

print("\nWORKER_2 Unmethylated List:", "\t", len(unmeth_list_worker_2), "\t", unmeth_list_worker_2[0])
print("\nWORKER_2 Methylated List:", "\t", len(meth_list_worker_2), "\t", meth_list_worker_2[0])

unmeth_list_drone = CreateDataSet("/Users/kristinaulicna/Documents/Rotation_1/Data/CpG_all_Drone_meth_less_0_1.txt", 1, 500)
meth_list_drone = CreateDataSet("/Users/kristinaulicna/Documents/Rotation_1/Data/CpG_all_Drone_meth_over_0_9.txt", 1, 319)

print("\nDRONE Unmethylated List:", "\t", len(unmeth_list_drone), "\t", unmeth_list_drone[0])
print("\nDRONE Methylated List:", "\t", len(meth_list_drone), "\t", meth_list_drone[0])


#---------- Create single list with 5 lists, each with 500 or 272 positions for [0], [1], [2], [3], [4] value:

from itertools import zip_longest
from FileHandling_Functions import BoxplotPlottingData
unmeth_pos_worker_1, meth_pos_worker_1 = BoxplotPlottingData(unmeth_list_worker_1, meth_list_worker_1)
unmeth_pos_worker_2, meth_pos_worker_2 = BoxplotPlottingData(unmeth_list_worker_2, meth_list_worker_2)
unmeth_pos_drone, meth_pos_drone = BoxplotPlottingData(unmeth_list_drone, meth_list_drone)


#---------- Figure out the normalisation factor (for both WORKER-2 and DRONE:

import numpy as np

index = 0
mean_list_drone = []
while index < 5:
    mean_list_drone.append(np.mean(unmeth_pos_drone[index]))
    index += 1
norm_factor_drone = abs(np.mean(mean_list_drone))
print ("\nDRONE Normalisation factor:", norm_factor_drone, type(norm_factor_drone), "\n")

index = 0
mean_list_work2 = []
while index < 5:
    mean_list_work2.append(np.mean(unmeth_pos_worker_2[index]))
    index += 1
norm_factor_work2 = abs(np.mean(mean_list_work2))
print ("\nWORK2 Normalisation factor:", norm_factor_work2, type(norm_factor_work2), "\n")


#---------- Normalise the WORKER-2 and DRONE data:

from FileHandling_Functions import NormaliseData
unmeth_pos_drone = NormaliseData(unmeth_pos_drone, norm_factor_drone)
meth_pos_drone = NormaliseData(meth_pos_drone, norm_factor_drone)

unmeth_pos_worker_2 = NormaliseData(unmeth_pos_worker_2, norm_factor_work2)
meth_pos_worker_2 = NormaliseData(meth_pos_worker_2, norm_factor_work2)


#---------- Zip the data according to the longest:
"""
def PrintZippedData(unmeth_pos, meth_pos):
    for un_pos, me_pos in zip_longest(unmeth_pos, meth_pos, fillvalue=None):
        index = 0
        while index < 5:
            print("WORKER_1 Unmeth_Pos:\t", len(unmeth_pos[index]), unmeth_pos[index][0], end="\t")
            print("WORKER_1 Meth_Pos:\t", len(meth_pos[index]), meth_pos[index][0])
            index += 1
        break

PrintZippedData(unmeth_pos_worker_1, meth_pos_worker_1)
PrintZippedData(unmeth_pos_worker_2, meth_pos_worker_2)
PrintZippedData(unmeth_pos_drone, meth_pos_drone)
"""





#----------- Plot a boxplot:

import matplotlib.pyplot as plt
import numpy as np

ticks = ["| _ _ _ _ C |", "| _ _ _ C G |", "| _ _ C G _ |", "| _ C G _ _ |", "| C G _ _ _ |"]

def set_box_color(bp, color):
    plt.setp(bp['boxes'], color=color)
    plt.setp(bp['whiskers'], color=color)
    plt.setp(bp['caps'], color=color)
    plt.setp(bp['medians'], color=color)

plt.figure(figsize=(6, 12))

    # Plot no.1:
plt.subplot(3, 1, 1)
bp_1_U = plt.boxplot(unmeth_pos_worker_1, positions=np.array(range(len(unmeth_pos_worker_1))) * 2.0 - 0.25, sym='', widths=0.4, patch_artist=True, boxprops=dict(alpha=.5))
bp_1_M = plt.boxplot(meth_pos_worker_1, positions=np.array(range(len(unmeth_pos_worker_1))) * 2.0 + 0.25, sym='', widths=0.4, patch_artist=True, boxprops=dict(alpha=.5))
plt.title('A) Worker 1')
plt.grid(b=None, axis='y')
plt.ylim(-25, 25)
plt.xlim(-1, len(ticks)*2 - 1)
plt.xticks(range(0, len(ticks) * 2, 2), ticks)
plt.plot([], c='coral', label='Unmethylated CpGs', linewidth = 3.0)
plt.plot([], c='c', label='Methylated CpGs', linewidth = 3.0)
plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=2)

    # Plot no.2:
plt.subplot(3, 1, 2)
bp_2_U = plt.boxplot(unmeth_pos_worker_2, positions=np.array(range(len(unmeth_pos_worker_1))) * 2.0 - 0.25, sym='', widths=0.4, patch_artist=True, boxprops=dict(alpha=.5))
bp_2_M = plt.boxplot(meth_pos_worker_2, positions=np.array(range(len(unmeth_pos_worker_1))) * 2.0 + 0.25, sym='', widths=0.4, patch_artist=True, boxprops=dict(alpha=.5))
plt.title('B) Worker 2')
plt.ylabel('Observed vs. Expected signal [pAmp] from Nanopush')
plt.grid(b=None, axis='y')
plt.ylim(-25, 25)
plt.xlim(-1, len(ticks)*2 - 1)
plt.xticks(range(0, len(ticks) * 2, 2), ticks)
plt.plot([], c='dodgerblue', label='Unmethylated CpGs', linewidth = 3.0)
plt.plot([], c='darkgoldenrod', label='Methylated CpGs', linewidth = 3.0)
plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=2)

    # Plot no.3:
plt.subplot(3, 1, 3)
bp_3_U = plt.boxplot(unmeth_pos_drone, positions=np.array(range(len(unmeth_pos_worker_1))) * 2.0 - 0.25, sym='', widths=0.4, patch_artist=True, boxprops=dict(alpha=.5))
bp_3_M = plt.boxplot(meth_pos_drone, positions=np.array(range(len(unmeth_pos_worker_1))) * 2.0 + 0.25, sym='', widths=0.4, patch_artist=True, boxprops=dict(alpha=.5))
plt.title('C) Drone')
plt.xlabel('Schematic of CpG in the pentamer', labelpad = 10)
plt.grid(b=None, axis='y')
plt.ylim(-25, 25)
plt.xlim(-1, len(ticks)*2 - 1)
plt.xticks(range(0, len(ticks) * 2, 2), ticks)
plt.plot([], c='orchid', label='Unmethylated CpGs', linewidth = 3.0)
plt.plot([], c='limegreen', label='Methylated CpGs', linewidth = 3.0)
plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=2)

    # Colours:
set_box_color(bp_1_U, 'coral')
set_box_color(bp_1_M, 'c')
set_box_color(bp_2_U, 'dodgerblue')
set_box_color(bp_2_M, 'darkgoldenrod')
set_box_color(bp_3_U, 'orchid')
set_box_color(bp_3_M, 'limegreen')

    # View:
plt.tight_layout()
plt.savefig("/Users/kristinaulicna/Documents/Rotation_1/Plots_For_Report/Boxplot_UandM.png")
plt.show()
plt.close()