def PlotSVMGamma(x_axis_data, y_axis_data, y_from, y_to, section=''):
    import matplotlib.pyplot as plt
    plt.figure()
    plt.semilogx(x_axis_data, y_axis_data)
    plt.title('semilogx')
    plt.ylim(y_from, y_to)
    plt.grid(b=None, axis='y')
    plt.tight_layout()
    plt.xlabel("Range of 'gamma' parameter")
    plt.ylabel("Sum of support vectors used for model training")
    plt.title("SVM Optimisation to the lowest number of SV used per category")
    plt.savefig("/Users/kristinaulicna/Documents/Rotation_1/Plots_For_Report/GammaSVM_" + section + ".png", bbox_inches = 'tight')
    plt.show()


#TODO: W O R K E R  1 = Training       . . . on 250 lines of RAW data (unmeth & meth)

    # Create data set by concatenating & reshuffling the unmeth & meth lists with their CpG methylation score:

import numpy as np
from FileHandling_Functions import CreateDataSet
from FileHandling_Functions import RandomiseLists

unmeth_train = CreateDataSet("/Users/kristinaulicna/Documents/Rotation_1/Data/CpG_all_Worker_1_meth_less_0_1.txt", 1, 250)  # out of 500 total
meth_train = CreateDataSet("/Users/kristinaulicna/Documents/Rotation_1/Data/CpG_all_Worker_1_meth_over_0_9.txt", 1, 250)    # out of 272 total
train_list, train_score = RandomiseLists(unmeth_train, meth_train)

unmeth_valid = CreateDataSet("/Users/kristinaulicna/Documents/Rotation_1/Data/CpG_all_Worker_1_meth_less_0_1.txt", 251, 272)
meth_valid = CreateDataSet("/Users/kristinaulicna/Documents/Rotation_1/Data/CpG_all_Worker_1_meth_over_0_9.txt", 251, 272)
valid_list, valid_score = RandomiseLists(unmeth_valid, meth_valid)


    # Train the SVM on Worker_1 CpGs (gamma = 10^-7 to 10^7):

from sklearn import svm
clf = svm.SVC(gamma = 'scale', kernel = 'rbf')
print (clf.fit(train_list, train_score))
print ("\nGamma = 'scale' : SV per class:", clf.n_support_)


    # Try the SVM with a range of 'gamma' parameters (10^-7 to 10^7):

parameter_range = 7
gamma_range = [10 ** i for i in range(-parameter_range, parameter_range)]

sv_sum_list = []
for gamma in gamma_range:
    clf = svm.SVC(gamma = gamma, kernel = 'rbf')
    clf.fit(train_list, train_score)
    sc = clf.score(train_list, train_score)
    sv = list(clf.n_support_)
    sv_sum = int(sv[0]) + int(sv[1])
    sv_sum_list.append(sv_sum)
    #print("Support vectors:", sv, "\t=", sv_sum, "for gamma =", gamma, "\t", "; score =", sc)


minimum = [i for i,x in enumerate(sv_sum_list) if x == min(sv_sum_list)]
print ("\nThe minimal sum of SV used is", min(sv_sum_list), "when 'gamma' is equal to", gamma_range[minimum[0]])
#PlotSVMGamma(gamma_range, sv_sum_list, 420, 520, section = 'full')


    # TODO: Plot Full Graph:
import matplotlib.pyplot as plt
plt.figure()

    # Main Figure:
plt.semilogx(gamma_range, sv_sum_list, c = 'k', color = 'coral', linewidth = 5.0)
plt.ylim(420, 520)
plt.grid(b=None, axis='y')
plt.tight_layout()
plt.xlabel("Range of 'gamma' parameter")
plt.ylabel("Sum of support vectors used for model training")

    # Fake the Legend:
plt.plot([], c='coral', linewidth = 5.0, label='SV number (full range)')
plt.plot([], c='c', linewidth = 2.0, label='SV number (selection)')
plt.legend(loc='upper center', ncol=2)


    # Try the SVM with a smaller range of 'gamma' parameters (0.0001 - 0.01):

gamma_range = [gamma_range[minimum[0]] * 10 ** -1, gamma_range[minimum[0]], gamma_range[minimum[0]] * 10 ** 1]
gamma_range = [float(i) / 10000 for i in list(range(1, 101, 1))]

sv_sum_list = []
for gamma in gamma_range:
    clf = svm.SVC(gamma = gamma, kernel = 'rbf')
    clf.fit(train_list, train_score)
    sc = clf.score(train_list, train_score)
    sv = list(clf.n_support_)
    sv_sum = int(sv[0]) + int(sv[1])
    sv_sum_list.append(sv_sum)
    print ("Support vectors:", sv, "=", sv_sum, "for gamma =", gamma, "\t", "; score =", sc)

minimum = [i for i,x in enumerate(sv_sum_list) if x == min(sv_sum_list)]
print ("\nThe minimal sum of SV used is", min(sv_sum_list), "when 'gamma' is equal to", gamma_range[minimum[0]])
#PlotSVMGamma(gamma_range[9:(len(sv_sum_list)-1)], sv_sum_list[9:(len(sv_sum_list)-1)], 420, 440, section = 'zoom')


    # Detailed Figure:
sub_axes = plt.axes([0.63, 0.23, 0.3, 0.3])
sub_axes.semilogx(gamma_range[9:(len(sv_sum_list)-1)], sv_sum_list[9:(len(sv_sum_list)-1)], c = 'k', color = 'c', linewidth = 2.0)
sub_axes.grid(b=None, which='major', axis='y')

    # Save & View
plt.savefig("/Users/kristinaulicna/Documents/Rotation_1/Plots_For_Report/GammaSVM.png", bbox_inches='tight')
plt.show()
plt.close()


    # Re-train the SVM with new 'gamma':

clf = svm.SVC(kernel='rbf', gamma=gamma_range[minimum[0]], C=1).fit(train_list, train_score)
print(clf.fit(train_list, train_score))



    #TODO: W O R K E R  1 = Validation     . . . on 22 lines of RAW data (unmeth & meth)

# 1.) Call the 'cross_val_score' helper function:

from sklearn.model_selection import cross_val_score
scores = cross_val_score(clf, valid_list, valid_score, cv=5, scoring='f1_macro')
print ("\nAccuracy (method 1): %0.2f (+/- %0.2f)" % (scores.mean(), scores.std() * 2))


# 2.) Pass a cross validation iterator:

from sklearn.model_selection import ShuffleSplit
cv = ShuffleSplit(n_splits=5, test_size=1, random_state=0)
scores = cross_val_score(clf, train_list, train_score, cv=cv)
print ("\nAccuracy (method 2): %0.2f (+/- %0.2f)" % (scores.mean(), scores.std() * 2))


# 3.) Use an iterable yielding (train, test) splits as arrays of indices

def custom_cv_2folds(data):
    n = len(data)
    i = 1
    while i <= 2:
        idx = np.arange(n * (i - 1) / 2, n * i / 2, dtype=int)
        yield idx, idx
        i += 1

custom_cv = custom_cv_2folds(valid_list)
scores = cross_val_score(clf, valid_list, valid_score, cv=custom_cv)
print ("\nAccuracy (method 3): %0.2f (+/- %0.2f)" % (scores.mean(), scores.std() * 2))


    # Define an SVM function to be called:

def TrainSVMonWorker():
    clf = svm.SVC(kernel = 'rbf', gamma = gamma_range[minimum[0]], C = 1).fit(train_list, train_score)
    print ("\nTrained Classifier (Worker_1):", clf.fit(train_list, train_score))
    print ("Number of SVs per category:", clf.n_support_)
    return clf
