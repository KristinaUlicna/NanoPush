import matplotlib.pyplot as plt
fig = plt.figure()
fig.text(0.5, 0.0, 'Nanopush (Worker-1-trained SVM) methylation score', ha='center')
fig.text(0.0, 0.5, 'Reference (WGBSeq or Nanopolish) methylation score', va='center', rotation='vertical')


    # Worker-1:

axis_Ref = []
axis_SVM = []
for counter, line in enumerate(open("/Users/kristinaulicna/Documents/Rotation_1/Data/Results_Worker_1.txt", "r")):
    line = line.strip().split('\t')
    axis_Ref.append(float(line[5]))
    axis_SVM.append(float(line[8]))

plt.subplot(2, 2, 1)
plt.scatter(axis_SVM, axis_Ref, alpha=0.2, color='darksalmon')
plt.xlim(-0.05, 1.05)
plt.ylim(-0.05, 1.05)
plt.title('A) Worker 1')


    # Worker-2:

axis_Ref = []
axis_SVM = []
for counter, line in enumerate(open("/Users/kristinaulicna/Documents/Rotation_1/Data/Results_Worker_2.txt", "r")):
    line = line.strip().split('\t')
    axis_Ref.append(float(line[5]))
    axis_SVM.append(float(line[8]))

plt.subplot(2, 2, 2)
plt.scatter(axis_SVM, axis_Ref, alpha=0.2, color='skyblue')
plt.xlim(-0.05, 1.05)
plt.ylim(-0.05, 1.05)
plt.title('B) Worker 2')


    # Drone-Nanopolish:

axis_Ref = []
axis_SVM = []
for counter, line in enumerate(open("/Users/kristinaulicna/Documents/Rotation_1/Data/Results_Drone_Polish.txt", "r")):
    line = line.strip().split('\t')
    axis_Ref.append(float(line[5]))
    axis_SVM.append(float(line[8]))

plt.subplot(2, 2, 3)
plt.scatter(axis_SVM, axis_Ref, alpha=0.2, color='plum')
plt.xlim(-0.05, 1.05)
plt.ylim(-0.05, 1.05)
plt.title('C) Drone (Nanopolish)')


    # Drone-WGBSeq:

axis_Ref = []
axis_SVM = []
for counter, line in enumerate(open("/Users/kristinaulicna/Documents/Rotation_1/Data/Results_Drone_WGBSeq.txt", "r")):
    line = line.strip().split('\t')
    axis_Ref.append(float(line[5]))
    axis_SVM.append(float(line[8]))

plt.subplot(2, 2, 4)
plt.scatter(axis_SVM, axis_Ref, alpha=0.2, color='yellowgreen')
plt.xlim(-0.05, 1.05)
plt.ylim(-0.05, 1.05)
plt.title('D) Drone (WGB-Seq)')


    # Save & View:

plt.tight_layout()
plt.savefig("/Users/kristinaulicna/Documents/Rotation_1/Plots_For_Report/Scatter_RefVsSVM.png", bbox_inches = 'tight')
plt.show()
plt.close()
