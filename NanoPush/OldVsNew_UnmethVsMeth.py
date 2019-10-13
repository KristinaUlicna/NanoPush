# Plot the boxplots for 4000 unmethylated & methylated CpGs each, analysed by old (Nanopolish) & new (NN events) method
# with their 5 values and see if there is a difference:

# Organise your data in the style of 100:5 matrix
# (100 rows corresponding to individual CpG & 5 columns corresponding to 5 signal data points per CpG)

un_old = [[] for i in range(5)]
for CpG in open("/Users/kristinaulicna/Documents/Rotation_1/Archive/DRONE/UnmethOldNoNone.txt"):
    CpG = CpG.strip().split('\t')
    CpG = [float(i) for i in CpG[6:11]]
    index = 0
    for number in CpG:
        un_old[index].append(CpG[index])
        index += 1
print ("Un_Old:", type(un_old[0][0]), un_old[0][0], un_old[1][0], un_old[2][0], un_old[3][0], un_old[4][0])


un_new = [[] for i in range(5)]
for CpG in open("/Users/kristinaulicna/Documents/Rotation_1/Archive/DRONE/UnmethNewNoNone.txt"):
    CpG = CpG.strip().split('\t')
    CpG = [float(i) for i in CpG[6:11]]
    index = 0
    for number in CpG:
        un_new[index].append(CpG[index])
        index += 1
print ("Un_New:", type(un_new[0][0]), un_new[0][0], un_new[1][0], un_new[2][0], un_new[3][0], un_new[4][0])


me_old = [[] for i in range(5)]
for CpG in open("/Users/kristinaulicna/Documents/Rotation_1/Archive/DRONE/MethOldNoNone.txt"):
    CpG = CpG.strip().split('\t')
    CpG = [float(i) for i in CpG[6:11]]
    index = 0
    for number in CpG:
        me_old[index].append(CpG[index])
        index += 1
print ("Me_Old:", type(me_old[0][0]), me_old[0][0], me_old[1][0], me_old[2][0], me_old[3][0], me_old[4][0])


me_new = [[] for i in range(5)]
for CpG in open("/Users/kristinaulicna/Documents/Rotation_1/Archive/DRONE/MethNewNoNone.txt"):
    CpG = CpG.strip().split('\t')
    CpG = [float(i) for i in CpG[6:11]]
    index = 0
    for number in CpG:
        me_new[index].append(CpG[index])
        index += 1
print ("Me_New:", type(me_new[0][0]), me_new[0][0], me_new[1][0], me_new[2][0], me_new[3][0], me_new[4][0])

print ("\nDone!")


#--------- Prepare the data:

import numpy as np

boxplot_datasize = 1500
all_lists = [un_old, un_new, me_old, me_new]

for single_list in all_lists:
    index = 0
    while index <= 4:
        single_list[index] = single_list[index][0:boxplot_datasize]
        index += 1

un_old = all_lists[0]
un_new = all_lists[1]
me_old = all_lists[2]
me_new = all_lists[3]

print ()
print ("Un_Old:", type(un_old[0][0]), len(un_old[0]), len(un_old[1]), len(un_old[2]), len(un_old[3]), len(un_old[4]), un_old[0][0], un_old[1][0], un_old[2][0], un_old[3][0], un_old[4][0])
print ("Un_New:", type(un_new[0][0]), len(un_new[0]), len(un_new[1]), len(un_new[2]), len(un_new[3]), len(un_new[4]), un_new[0][0], un_new[1][0], un_new[2][0], un_new[3][0], un_new[4][0])
print ("Me_Old:", type(me_old[0][0]), len(me_old[0]), len(me_old[1]), len(me_old[2]), len(me_old[3]), len(me_old[4]), me_old[0][0], me_old[1][0], me_old[2][0], me_old[3][0], me_old[4][0])
print ("Me_New:", type(me_new[0][0]), len(me_new[0]), len(me_new[1]), len(me_new[2]), len(me_new[3]), len(me_new[4]), me_new[0][0], me_new[1][0], me_new[2][0], me_new[3][0], me_new[4][0])

mean_un_old = [np.mean(un_old[i]) for i in range(5)]
mean_me_old = [np.mean(me_old[i]) for i in range(5)]
mean_un_new = [np.mean(un_new[i]) for i in range(5)]
mean_me_new = [np.mean(me_new[i]) for i in range(5)]
print ("\nMean_Un_Old:", mean_un_old, "\t", "Mean_Me_Old:", mean_me_old, "\t", "Mean_Un_New:", mean_un_new, "\t", "Mean_Me_New:", mean_me_new)
print ("Mean_Un_New:", mean_un_new)
print ("Shift:", np.mean(mean_un_new))
    #TODO: Add 15+ to all data created with NEW method & train an SVM.

#un_new_shifted = [(number + 15) for number in un_new[index] for index in range(4)]
#me_new_shifted = [(number + 15) for number in me_new[index] for index in range(4)]
    #TODO: Ask Nunik how to make a list comprehension with 2 for loops.

un_new_shifted = [[] for i in range(5)]
me_new_shifted = [[] for i in range(5)]
for i in range(5):
    for number_un, number_me in zip(un_new[i], me_new[i]):
        un_new_shifted[i].append(number_un+15)
        me_new_shifted[i].append(number_me+15)

print ("\nUn_New + shift:", type(un_new_shifted[0][0]), len(un_new_shifted[0]), len(un_new_shifted[1]), len(un_new_shifted[2]), len(un_new_shifted[3]), len(un_new_shifted[4]), un_new_shifted[0][0], un_new_shifted[1][0], un_new_shifted[2][0], un_new_shifted[3][0], un_new_shifted[4][0])
print ("Me_New + shift:", type(me_new_shifted[0][0]), len(me_new_shifted[0]), len(me_new_shifted[1]), len(me_new_shifted[2]), len(me_new_shifted[3]), len(me_new_shifted[4]), me_new_shifted[0][0], me_new_shifted[1][0], me_new_shifted[2][0], me_new_shifted[3][0], me_new_shifted[4][0])


#---------- Create a file with normalised data (which have +15 added to them):

un_new_shifted_data_file = open("/Users/kristinaulicna/Documents/Rotation_1/Archive/DRONE/UnmethNewDataOnly.txt", "w")
string = ''
order = 0
while order < len(un_new_shifted[0]):
    string += (str(un_new_shifted[0][order]) + "\t" + str(un_new_shifted[1][order]) + "\t" + str(un_new_shifted[2][order]) + "\t" + str(un_new_shifted[3][order]) + "\t" + str(un_new_shifted[4][order]))
    string += "\n"
    order += 1
un_new_shifted_data_file.write(string)
un_new_shifted_data_file.close()

me_new_shifted_data_file = open("/Users/kristinaulicna/Documents/Rotation_1/Archive/DRONE/MethNewDataOnly.txt", "w")
string = ''
order = 0
while order < len(me_new_shifted[0]):
    string += (str(me_new_shifted[0][order]) + "\t" + str(me_new_shifted[1][order]) + "\t" + str(me_new_shifted[2][order]) + "\t" + str(me_new_shifted[3][order]) + "\t" + str(me_new_shifted[4][order]))
    string += "\n"
    order += 1
me_new_shifted_data_file.write(string)
me_new_shifted_data_file.close()


#---------- Plot a boxplot:

import matplotlib.pyplot as plt

ticks = ["| _ _ _ _ C |", "| _ _ _ C G |", "| _ _ C G _ |", "| _ C G _ _ |", "| C G _ _ _ |"]

def set_box_color(bp, color):
    plt.setp(bp['boxes'], color=color)
    plt.setp(bp['whiskers'], color=color)
    plt.setp(bp['caps'], color=color)
    plt.setp(bp['medians'], color=color)

plt.figure()
#TODO: Adjust the size of the figure to fit the title and axes labels!

bpl = plt.boxplot(un_old, positions=np.array(range(len(un_old)))*2.0-0.75, sym='', widths=0.4)
bpr = plt.boxplot(un_new, positions=np.array(range(len(un_new)))*2.0+0.25, sym='', widths=0.4)
bpll = plt.boxplot(me_old, positions=np.array(range(len(me_old)))*2.0-0.25, sym='', widths=0.4)
bprr = plt.boxplot(me_new, positions=np.array(range(len(me_new)))*2.0+0.75, sym='', widths=0.4)
    # optional: set extra parameter as 'showmeans=True'

set_box_color(bpl, 'orange')
set_box_color(bpr, 'purple')
set_box_color(bpll, 'blue')
set_box_color(bprr, 'green')

# draw temporary lines and use them to create a legend
plt.plot([], c='orange', label='Unmethylated OLD')
plt.plot([], c='blue', label='Methylated OLD')
plt.plot([], c='purple', label='Unmethylated NEW')
plt.plot([], c='green', label='Methylated NEW')
plt.legend()


plt.xticks(range(0, (len(ticks) * 2) + 2, 2), ticks)
plt.xlim(-2, (len(ticks) * 2))
plt.ylim(-70, 30)
plt.grid(b=None, axis='y')
plt.tight_layout()
plt.xlabel("CpG schematic in the pentamer")
plt.ylabel("Observed vs. Expected signal value")
plt.title("D R O N E : "
          "\nComparison of UNMETH (orange, purple) & METH (blue, green)"
          "\nCpG signals using OLD (orange, blue) & NEW (purple, green)"
          "\nmethod across the whole pentamer")
plt.savefig("/Users/kristinaulicna/Documents/Rotation_1/Archive/BoxplotComparison1500CpGs.png")
plt.show()
