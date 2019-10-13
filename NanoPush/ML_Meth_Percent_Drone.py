#TODO: Use SVM trained on WORKER data to predict the meth score of the DRONE.
#TODO: Append the list of drone data...
#TODO: Do those compare to the Nanopolish predictions?

from ML_SVM_Worker import TrainSVMonWorker
clf = TrainSVMonWorker()
print (clf)


    # Function to create a file with Rob's Nanopolish & my Nanopush meth score:

def ReferenceVsMySVM(iterable_file, reference_file, reference_file_combo):
    output_file = open(reference_file_combo, "w")
    input_file = open(reference_file, "r")
    CpG = list(input_file.readline().strip().split('\t')) + ['0', '0', '0']
    for line in open(iterable_file, "r"):
        line = line.strip().split('\t')
        data = [(float(i) + 15.97) for i in line[6:11]]
        #TODO: This normalisation factor does not apply to both Worker-2 and Drone!
        if line[2] == '-':
            line[1] = int(line[1]) - 1
        if line[0] == CpG[0]:
            if line[0] == CpG[0] and int(line[1]) == (int(CpG[1])-1):
                CpG[6] = int(CpG[6]) + 1
                if int(clf.predict([data])) == 1:
                    CpG[7] = int(CpG[7]) + 1
            elif line[0] == CpG[0] and int(line[1]) > (int(CpG[1])-1):
                if int(CpG[6]) == 0:
                    CpG[8] = float(CpG[8])
                else:
                    CpG[8] = int(CpG[7]) / int(CpG[6])
                #print(CpG)
                if int(CpG[6]) >= 5:
                    string = ''
                    for parameter in CpG:
                        string += (str(parameter) + "\t")
                    string = string[:-1]
                    string += "\n"
                    output_file.write(string)
                CpG = list(input_file.readline().strip().split('\t')) + ['0', '0', '0']
                if line[0] == CpG[0] and int(line[1]) == (int(CpG[1])-1):
                    CpG[6] = int(CpG[6]) + 1
                    if int(clf.predict([data])) == 1:
                        CpG[7] = int(CpG[7]) + 1
            elif line[0] == CpG[0] and int(line[1]) < (int(CpG[1])-1):
                continue
        else:
            CpG = list(input_file.readline().strip().split('\t')) + ['0', '0', '0']


ReferenceVsMySVM("/Users/kristinaulicna/Documents/Rotation_1/Data/Sorted_Worker_1_NoNone.txt",\
                    "/Users/kristinaulicna/Documents/Rotation_1/Data/Reference_WGBSeq.txt",\
                    "/Users/kristinaulicna/Documents/Rotation_1/Data/Results_Worker_1.txt")

ReferenceVsMySVM("/Users/kristinaulicna/Documents/Rotation_1/Data/Sorted_Worker_2_NoNone.txt",\
                    "/Users/kristinaulicna/Documents/Rotation_1/Data/Reference_WGBSeq.txt",\
                    "/Users/kristinaulicna/Documents/Rotation_1/Data/Results_Worker_2.txt")

ReferenceVsMySVM("/Users/kristinaulicna/Documents/Rotation_1/Data/Sorted_Drone_NoNone.txt",\
                    "/Users/kristinaulicna/Documents/Rotation_1/Data/Reference_Nanopolish.txt",\
                    "/Users/kristinaulicna/Documents/Rotation_1/Data/Results_Drone_Polish.txt")

ReferenceVsMySVM("/Users/kristinaulicna/Documents/Rotation_1/Data/Sorted_Drone_NoNone.txt",\
                    "/Users/kristinaulicna/Documents/Rotation_1/Data/Reference_WGBSeq.txt",\
                    "/Users/kristinaulicna/Documents/Rotation_1/Data/Results_Drone_WGBSeq.txt")


# Plot the Scatter plot:

def PlotScatterSVMvsRef(meth_percentage_file, which_dataset):
    axis_Ref = []
    axis_SVM = []
    for counter, line in enumerate(open(meth_percentage_file, "r")):
        line = line.strip().split('\t')
        axis_Ref.append(float(line[5]))
        axis_SVM.append(float(line[8]))
    import matplotlib.pyplot as plt
    plt.scatter(axis_SVM, axis_Ref, alpha=0.2)
    plt.xlim(-0.05, 1.05)
    plt.ylim(-0.05, 1.05)
    plt.xlabel("My Worker-1-trained SVM on " + which_dataset)
    plt.ylabel("Rob's Reference (WGBSeq or Nanopolish)" + which_dataset)
    plt.title("Difference in CpG methylation scoring\nwith Rob's reference & my SVM; " + which_dataset)
    plt.savefig("/Users/kristinaulicna/Documents/Rotation_1/Archive/Scatter_RefVsSVM_" + which_dataset + ".png", bbox_inches = 'tight')
    plt.show()

PlotScatterSVMvsRef("/Users/kristinaulicna/Documents/Rotation_1/Data/Results_Worker_1.txt", "WORKER_1")
PlotScatterSVMvsRef("/Users/kristinaulicna/Documents/Rotation_1/Data/Results_Worker_2.txt", "WORKER_2")
PlotScatterSVMvsRef("/Users/kristinaulicna/Documents/Rotation_1/Data/Results_Drone_Polish.txt", "DRONE_Polish")
PlotScatterSVMvsRef("/Users/kristinaulicna/Documents/Rotation_1/Data/Results_Drone_WGBSeq.txt", "DRONE_WGBSeq")


#TODO: May need to add +1 to the mitochondrial DNA!
#TODO: If filtered for the coverage > 5, the graph suggests there is not a single CpG scored over 0.9 in the Rob's method... REALLY?
