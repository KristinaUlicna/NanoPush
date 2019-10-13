#TODO: Make sure the interpreter for this script is set to python2.7!

import numpy as np
import keras
from keras.models import Sequential
from keras.layers import Dense, Dropout, Activation
from keras.optimizers import SGD
from keras.utils import to_categorical

data_list = []
for line in open("/Users/kristinaulicna/Documents/Rotation_1/CpG_Analysis/NN_training_dataset_MethScore.txt", "r"):
    line_list = line.strip().split('\t')
    data_list.append(line_list)
#print ("Data List:", len(data_list), data_list[0], data_list[1], type(data_list[1][0]))


# Nunik's code:
data = data_list[1:]
labels = np.array([data[i][-1] for i in range(len(data))])
features = np.array([data[i][:-1] for i in range(len(data))])
labels_categorical = to_categorical(labels)

print(labels.shape)
print(features.shape)

model = Sequential()
# Dense(64) is a fully-connected layer with 64 hidden units.
# in the first layer, you must specify the expected input data shape:
# here, 20-dimensional vectors.
model.add(Dense(64, activation='relu', input_dim=5))
model.add(Dropout(0.5))
model.add(Dense(64, activation='relu'))
model.add(Dropout(0.5))
model.add(Dense(2, activation='softmax'))

print(model.summary())

sgd = SGD(lr=0.01, decay=1e-6, momentum=0.9, nesterov=True)
model.compile(loss='categorical_crossentropy',
              optimizer=sgd,
              metrics=['accuracy'])

model.fit(features, labels_categorical,
          epochs=20,
          batch_size=8, validation_split=0.2)

score = model.evaluate(features, labels_categorical, batch_size=128)
print(score)
print("Prediction {}".format(model.predict(features[:1,:])))
print(labels[0])
print(labels_categorical[0,:])


# Nunka's code:
        # Now you can see if the NN gets the meth_score right: Pick 5 random lists for unmethylated CpGs...
        # TODO: WATCH OUT! This is python2.7!!! (That's why altered printing...)

try_1 = np.array([18.944628359199996, 54.876043791, 4.984124630099998, 12.705428988999998, 49.472132714])
try_1 = try_1.reshape((1, 5))
print (model.predict(try_1))
    # Expected: [[1. 0.]]; Observed: [[1. 0.]]
try_2 = np.array([3.061183359200001, 35.422293791, 13.206959369899991, 7.363685677666666, 26.43067471399999])
try_2 = try_2.reshape((1, 5))
print (model.predict(try_2))
    # Expected: [[1. 0.]]; Observed: [[1. 0.]]
try_3 = np.array([3.313106640800001, 40.205907124333336, 2.673138869900015, 7.585100988999997, 26.796412714])
try_3 = try_3.reshape((1, 5))
print (model.predict(try_3))
    # Expected: [[1. 0.]]; Observed: [[1. 0.]]
try_4 = np.array([9.017468359199995, 38.069533791, 4.7444940699, 2.6389459889999927, 31.638082714000006])
try_4 = try_4.reshape((1, 5))
print (model.predict(try_4))
    # Expected: [[1. 0.]]; Observed: [[1. 0.]]
try_5 = np.array([49.472132714, 17.816624521699993, 33.56380045326667, 8.3014805368, 27.99538970056667])
try_5 = try_5.reshape((1, 5))
print (model.predict(try_5))
   # Expected: [[1. 0.]]; Observed: [[1. 0.]]


# TODO: Convert one hot encoding into integers!

    # input:
a = np.array([[0,1,0,0],[1,0,0,0],[0,0,0,1]])
a = [np.where(r == 1)[0][0] for r in a]
print (a)
print ("END!")
    # output:   [1, 0, 3]
                # = in the 1st list [0], the '1' is on position '2' [1]...
                # = in the 2nd list [1], the '1' is on position '1' [0]...
                # = in the 3rd list [2], the '1' is on position '4' [3]...


# Loop through methylated lists of data in the 'over_0_9' txt file:

results_list = []
for order, values in enumerate(open("/Users/kristinaulicna/Documents/Rotation_1/CpG_Analysis/CpGs_Drone_meth_over_0_9_ML_vectors.txt", "r")):
    if order <= 100:
        continue
    else:
        values = values.strip().split('\t')
        values = np.array(values[0:5]).reshape((1, 5))
        results_list.append(list(model.predict(values)[0]))
        print (model.predict(values)[0])

results_list = np.array(results_list)
results_list = [np.where(r == 1)[0][0] for r in results_list]
print (len(results_list), results_list)

counter_right = 0
counter_wrong = 0

for element in results_list:
    if element == 0:
        counter_wrong += 1
    if element == 1:
        counter_right += 1
print ("Counter Right:", counter_right)
print ("Counter Wrong:", counter_wrong)