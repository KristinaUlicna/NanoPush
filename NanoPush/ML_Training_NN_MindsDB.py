# Start with Machine Learning!
        # Method:  N E U R A L  N E T W O R K S


# Prepare a file for training: 100 unmethylated & 100 methylated CpG values with their score in the 6th column, with a header!

def PrepareNNDataset(input_UNMETH_file, input_METH_file, output_file, how_long):
    NN_training_data = open(output_file, "w")
    header = ['CpG_position_5', 'CpG_position_4', 'CpG_position_3', 'CpG_position_2', 'CpG_position_1', 'meth_score']
    file_header = ''
    for word in header:
        file_header += word + "\t"
    file_header = file_header[:-1]
    file_header += "\n"
    NN_training_data.write(file_header)
    for line_order, line in enumerate(open(input_UNMETH_file, "r")):
        new_line = line.strip().split('\t')
        if len(new_line) == 11:
            new_line = new_line[6:11]
        new_line.append('0')
        if line_order > int(how_long):
            break
        string = ''
        for item in new_line:
            string += item + "\t"
        string = string[:-1]
        string += "\n"
        NN_training_data.write(string)
    for line_order, line in enumerate(open(input_METH_file, "r")):
        new_line = line.strip().split('\t')
        if len(new_line) == 11:
            new_line = new_line[6:11]
        new_line.append('1')
        if line_order > int(how_long):
            break
        string = ''
        for item in new_line:
            string += item + "\t"
        string = string[:-1]
        string += "\n"
        NN_training_data.write(string)
    NN_training_data.close()


PrepareNNDataset("/Users/kristinaulicna/Documents/Rotation_1/Archive/DRONE/NewSelectedFromRob/CpGs_sel_Drone_meth_less_0_1.txt", "/Users/kristinaulicna/Documents/Rotation_1/Archive/DRONE/NewSelectedFromRob/CpGs_sel_Drone_meth_over_0_9.txt", "/Users/kristinaulicna/Documents/Rotation_1/Archive/NN_Training_Data/NN_DRONE_Training_Data.txt", 200)
PrepareNNDataset("/Users/kristinaulicna/Documents/Rotation_1/Archive/WORKER/CpG_all_Worker_meth_less_0_1.txt", "/Users/kristinaulicna/Documents/Rotation_1/Archive/WORKER/CpG_all_Worker_meth_over_0_9.txt", "/Users/kristinaulicna/Documents/Rotation_1/Archive/NN_Training_Data/NN_WORKER_Testing_Data.txt", 200)


# Training & using the neural network model:

from mindsdb import *

# first we initiate MindsDB
mdb = MindsDB()

# we tell mindsDB what we want to learn and from what data
mdb.learn(
    from_data = "/Users/kristinaulicna/Documents/Rotation_1/Archive/NN_Training_Data/NN_DRONE_Training_Data.txt", # the path to the file where we can learn from, (note: can be url)
    predict = 'meth_score', # the column we want to learn to predict given all the data in the file
    model_name = 'meth_prediction' # the name of this model
)

# use the model to make predictions
#mdb.predict(predict='meth_score', when={'CpG_position_5' : 18.944628359199996, 'CpG_position_4' : 54.876043791, 'CpG_position_3' : 4.984124630099998, 'CpG_position_2' : 12.705428988999998, 'CpG_position_1' : 49.472132714}, model_name='meth_prediction')

#TODO: Loop through the data:
for order, line in enumerate(open("/Users/kristinaulicna/Documents/Rotation_1/Archive/NN_Training_Data/NN_WORKER_Testing_Data.txt", "r")):
    line = line.strip().split('\t')
    if order == 0:
        continue
    print ("Line:", order, "\t", mdb.predict(predict='meth_score', when={'CpG_position_5' : line[0], 'CpG_position_4' : line[1], 'CpG_position_3' : line[2], 'CpG_position_2' : line[3], 'CpG_position_1' : line[4]}, model_name='meth_prediction'))


# you can now print the results
#print('The predicted price is ${score} with {conf} confidence'.format(score=result.predicted_values[0]['meth_prediction'], conf=result.predicted_values[0]['prediction_confidence']))