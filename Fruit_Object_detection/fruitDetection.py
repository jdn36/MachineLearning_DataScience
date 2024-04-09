#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


@author: josue_nataren
"""


import matplotlib.pyplot as plt
import numpy as np
import os, os.path
import tensorflow as tf
from keras.constraints import max_norm
from keras.utils import plot_model
from PIL import Image
import xml.etree.ElementTree as ET
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, classification_report

def custom_loss(y_true, y_pred):
    # Splitting true and predicted values for classification and regression
    y_true_cls, y_true_reg = tf.split(y_true, [1, 4], axis=-1)
    y_pred_cls, y_pred_reg = tf.split(y_pred, [1, 4], axis=-1)

    # Classification loss (softmax cross-entropy)
    cls_loss = tf.keras.losses.sparse_categorical_crossentropy(y_true_cls, y_pred_cls)

    # Regression loss (mean squared error)
    reg_loss = tf.keras.losses.mean_squared_error(y_true_reg, y_pred_reg)

    # Combine the losses (you can adjust the weights as needed)
    total_loss = cls_loss + reg_loss

    return total_loss

# Build your model
model = tf.keras.models.Sequential([
    # Define your model layers
])

def resize_image(input_path, new_width, new_height):
    # Open the image file
    image = Image.open(input_path)
    
    original_size = image.size
    
    resized_image = image.resize((new_width, new_height))

    # Convert the resized image to RGB mode (three channels)
    resized_image_rgb = resized_image.convert('RGB')
    
    # Convert the RGB image to a NumPy array
    resized_array = np.array(resized_image_rgb)
    
    
    return resized_array, original_size
    

# Step 1: Parse XML Files
def parse_xml(xml_file):
    tree = ET.parse(xml_file)
    root = tree.getroot()
    labels = []
    bboxes = []
    for obj in root.findall('object'):
        label = obj.find('name').text
        bbox = obj.find('bndbox')
        xmin = int(bbox.find('xmin').text)
        ymin = int(bbox.find('ymin').text)
        xmax = int(bbox.find('xmax').text)
        ymax = int(bbox.find('ymax').text)
        labels.append(label)
        bboxes.append((xmin, ymin, xmax, ymax))
    return labels, bboxes

train_folder_path = 'Fruits_OD/train_zip/train/'
test_folder_path = 'Fruits_OD/test_zip/test/'


ClassificationLabel = ["apple","banana","orange","mixed"]
#ClassificationLabel = [0,1,2,3]

#This is for A cases
totalFiles = 0

names = []
count = 0

allfilesA = []
outputLayer = list()

word_to_find = "mixed"

sizeOfImg = 300

vectorImages = np.zeros((0, sizeOfImg, sizeOfImg, 3))
outputLayer = np.zeros((0, 7))

dummyCounter = 0
currentX = 0
currentY = 0

for base, dirs, files in os.walk(train_folder_path):
    # Iterate over the files
    #print("We are here: ",type(files))
    files.sort()
    
    for file in files:
        # Check if the file has a JPEG or JPG extension
        #print(file,'\n')
        
        if word_to_find.lower() in file.lower():
            print("not using files: ",word_to_find)
        else:
            if file.lower().endswith(('.jpeg', '.jpg')):
                # If yes, construct the full file path and add it to the list
                file_path = os.path.join(base, file)
                allfilesA.append(file_path)
                dslocal, original_size = resize_image(file_path, sizeOfImg, sizeOfImg)
                # print(dslocal.shape)
                vectorImages = np.append(vectorImages, [dslocal], axis=0)
                # print(list(dslocal.shape))
                currentX = original_size[0]
                currentY = original_size[1]
                
                print(file_path)
            elif file.lower().endswith(('.xml')):
                labels, bboxes = parse_xml(train_folder_path+file)
                # print(labels)
                # print(bboxes)
                
                bndWidth = int( ( ( list(bboxes[0])[2] - list(bboxes[0])[0] )/currentX )*sizeOfImg )
                bndHeight = int( ( ( list(bboxes[0])[3] - list(bboxes[0])[1] )/currentY )*sizeOfImg )
                xCenter = int( ( ( list(bboxes[0])[2] + list(bboxes[0])[0] )/ (currentX*2) )*sizeOfImg )
                yCenter = int( ( ( list(bboxes[0])[3] + list(bboxes[0])[1] )/ (currentY*2) )*sizeOfImg )
                
                
                # bndWidth = int( ( ( list(bboxes[0])[0] )/currentX )*sizeOfImg )
                # bndHeight = int( ( ( list(bboxes[0])[1] )/currentY )*sizeOfImg )
                # xCenter = int( ( ( list(bboxes[0])[2] )/ (currentX) )*sizeOfImg )
                # yCenter = int( ( ( list(bboxes[0])[3] )/ (currentY) )*sizeOfImg )
                
                
                
                if labels[0] == "apple":
                    localArray = np.array([1, 0, 0, xCenter, yCenter, bndWidth, bndHeight])
                    print(localArray)
                    outputLayer = np.append(outputLayer, [localArray], axis=0)
                elif labels[0] == "banana":
                    localArray = np.array([0, 1, 0, xCenter, yCenter, bndWidth, bndHeight])
                    print(localArray)
                    outputLayer = np.append(outputLayer, [localArray], axis=0)
                    # outputLayer.append([1, 1, list(bboxes[0])])
                elif labels[0] == "orange":
                    localArray = np.array([0, 0, 1, xCenter, yCenter, bndWidth, bndHeight])
                    print(localArray)
                    outputLayer = np.append(outputLayer, [localArray], axis=0)
                    # outputLayer.append([1, 2, list(bboxes[0])])
       
            
       
        
       
        
       
        
       






#IF ADDING TEST TO THE DATA JUST TO CHECK IF IT IMPROVES

for base, dirs, files in os.walk(test_folder_path):
    # Iterate over the files
    #print("We are here: ",type(files))
    files.sort()
    
    for file in files:
        # Check if the file has a JPEG or JPG extension
        #print(file,'\n')
        
        if word_to_find.lower() in file.lower():
            print("not using files: ",word_to_find)
        else:
            if file.lower().endswith(('.jpeg', '.jpg')):
                # If yes, construct the full file path and add it to the list
                file_path = os.path.join(base, file)
                dslocal, original_size = resize_image(file_path, sizeOfImg, sizeOfImg)
                # print(dslocal.shape)
                vectorImages = np.append(vectorImages, [dslocal], axis=0)
                # print(list(dslocal.shape))
                currentX = original_size[0]
                currentY = original_size[1]
                
                print(file_path)
            elif file.lower().endswith(('.xml')):
                labels, bboxes = parse_xml(test_folder_path+file)
                # print(labels)
                # print(bboxes)
                
                
                bndWidth = int( ( ( list(bboxes[0])[2] - list(bboxes[0])[0] )/currentX )*sizeOfImg )
                bndHeight = int( ( ( list(bboxes[0])[3] - list(bboxes[0])[1] )/currentY )*sizeOfImg )
                xCenter = int( ( ( list(bboxes[0])[2] + list(bboxes[0])[0] )/ (currentX*2) )*sizeOfImg )
                yCenter = int( ( ( list(bboxes[0])[3] + list(bboxes[0])[1] )/ (currentY*2) )*sizeOfImg )
                
                
                # bndWidth = int( ( ( list(bboxes[0])[0] )/currentX )*sizeOfImg )
                # bndHeight = int( ( ( list(bboxes[0])[1] )/currentY )*sizeOfImg )
                # xCenter = int( ( ( list(bboxes[0])[2] )/ (currentX) )*sizeOfImg )
                # yCenter = int( ( ( list(bboxes[0])[3] )/ (currentY) )*sizeOfImg )
                
                
                if labels[0] == "apple":
                    localArray = np.array([1, 0, 0, xCenter, yCenter, bndWidth, bndHeight])
                    print(localArray)
                    outputLayer = np.append(outputLayer, [localArray], axis=0)
                elif labels[0] == "banana":
                    localArray = np.array([0, 1, 0, xCenter, yCenter, bndWidth, bndHeight])
                    print(localArray)
                    outputLayer = np.append(outputLayer, [localArray], axis=0)
                    # outputLayer.append([1, 1, list(bboxes[0])])
                elif labels[0] == "orange":
                    localArray = np.array([0, 0, 1, xCenter, yCenter, bndWidth, bndHeight])
                    print(localArray)
                    outputLayer = np.append(outputLayer, [localArray], axis=0)
                    # outputLayer.append([1, 2, list(bboxes[0])])
       







print("Training X shape: "+str(vectorImages.shape))
print("Training Y shape: "+str(outputLayer.shape))






x_data = vectorImages
y_data = outputLayer

x_train, x_validate, y_train, y_validate = train_test_split(x_data, y_data, test_size=0.2, random_state=42)



epochsToUse = 200
batchSizeToUse = 16

opt = tf.keras.optimizers.Adam(learning_rate=0.0005)





# input_layer = tf.keras.layers.Input((sizeOfImg, sizeOfImg, 3))

# # Convolutional layers
# conv1 = tf.keras.layers.Conv2D(32, kernel_size = (3,3), activation = 'tanh', padding='same')(input_layer)
# pool1 = tf.keras.layers.MaxPooling2D((2, 2))(conv1)
# conv2 = tf.keras.layers.Conv2D(64, kernel_size = (3,3), activation = 'tanh', padding='same')(pool1)
# pool2 = tf.keras.layers.MaxPooling2D((2, 2))(conv2)
# conv3 = tf.keras.layers.Conv2D(128, kernel_size = (3,3), activation = 'tanh', padding='same')(pool2)
# pool3 = tf.keras.layers.MaxPooling2D((2, 2))(conv3)
# # conv4 = tf.keras.layers.Conv2D(256, filters = 30, kernel_size = (3,3), activation = 'tanh', padding='same')(pool3)
# # pool4 = tf.keras.layers.MaxPooling2D((2, 2))(conv4)

# # Flatten layer
# flatten_layer = tf.keras.layers.Flatten()(pool3)

# # Shared dense layers
# # dense1 = tf.keras.layers.Dense(256, activation='relu')(flatten_layer)
# dense2 = tf.keras.layers.Dense(128, activation='relu')(flatten_layer)
# dense3 = tf.keras.layers.Dense(64, activation='relu')(dense2)
# dense4 = tf.keras.layers.Dense(32, activation='relu')(dense3)
# dense5 = tf.keras.layers.Dense(16, activation='relu')(dense4)

# # Classification output branch
# # output_classification = tf.keras.layers.Dense(1, activation=None, name='classification')(dense5)
# output_classification = tf.keras.layers.Dense(3, activation='softmax', name='classification')(dense5)

# # Regression output branch
# output_regression = tf.keras.layers.Dense(4, activation='linear', name='regression')(dense5)  # Assuming 4 values for bounding box coordinates

# # Define the model with multiple output branches
# model = tf.keras.models.Model(inputs=input_layer, outputs=[output_classification, output_regression])

# model.summary()


# # Compile the model
# model.compile(optimizer='adam',
#               loss={'classification': 'categorical_crossentropy',  # Use sparse_categorical_crossentropy since y_train contains class indices
#                     'regression': 'mse'},  # Use mean squared error for regression
#               metrics={'classification': 'accuracy', 'regression': 'mse'})  # Define metrics for each output branch

# # Fit the model
# history1 = model.fit(x_train, {'classification': y_train[:, 0:3], 'regression': y_train[:, 3:]},  # Splitting y_train for classification and regression outputs
#                     validation_data=(x_validate, {'classification': y_validate[:, 0:3], 'regression': y_validate[:, 3:]}),  # Splitting y_validate for validation
#                     epochs=epochsToUse,
#                     batch_size=batchSizeToUse)







# model = tf.keras.models.Sequential()

input_layer = tf.keras.layers.Input((sizeOfImg, sizeOfImg, 3))
conv1 = tf.keras.layers.Conv2D(kernel_size = (3,3), filters = 30,activation = 'tanh')(input_layer)
conv2 = tf.keras.layers.Conv2D(filters = 30, kernel_size = (3,3), activation = 'tanh')(conv1)
pool1 = tf.keras.layers.MaxPool2D(2,2)(conv2)
conv3 = tf.keras.layers.Conv2D(filters = 30, kernel_size = (3,3), activation = 'tanh')(pool1)
pool2 = tf.keras.layers.MaxPool2D(3,3)(conv3)
flatten1 = tf.keras.layers.Flatten()(pool2)
dense1 = tf.keras.layers.Dense(20, activation = 'relu')(flatten1)
dense2 = tf.keras.layers.Dense(15, activation = 'relu')(dense1)

output_classification = tf.keras.layers.Dense(3, activation='softmax', name='classification')(dense2)
output_regression = tf.keras.layers.Dense(4, activation='linear', name='regression')(dense2)  # Assuming 4 values for bounding box coordinates








# #ONLY CLASSIFICATION
# #Define the model with multiple output branches
# model = tf.keras.models.Model(inputs=input_layer, outputs=output_classification)
# # model = tf.keras.models.Model(inputs=input_layer, outputs=[output_classification, output_regression])

# model.summary()
# # Compile the model
# model.compile(optimizer=opt,
#               loss='categorical_crossentropy',  # Use mean squared error for regression
#               metrics='accuracy')  # Define metrics for each output branch
# # Fit the model
# history1 = model.fit(x_train, y_train[:, 0:3], 
#                     validation_data=(x_validate, y_validate[:, 0:3]),
#                     epochs=epochsToUse,
#                     batch_size=batchSizeToUse)






#CLASSIFICATION AND DETECTION
model = tf.keras.models.Model(inputs=input_layer, outputs=[output_classification, output_regression])

model.summary()
# Compile the model
model.compile(optimizer=opt,
              loss={'classification': 'categorical_crossentropy',  # Use sparse_categorical_crossentropy since y_train contains class indices
                    'regression': 'mse'},  # Use mean squared error for regression
              metrics={'classification': 'accuracy', 'regression': 'mse'})  # Define metrics for each output branch

# Fit the model
history1 = model.fit(x_train, {'classification': y_train[:, 0:3], 'regression': y_train[:, 3:]},  # Splitting y_train for classification and regression outputs
                    validation_data=(x_validate, {'classification': y_validate[:, 0:3], 'regression': y_validate[:, 3:]}),  # Splitting y_validate for validation
                    epochs=epochsToUse,
                    batch_size=batchSizeToUse)






# opt = tf.keras.optimizers.Adam(learning_rate=0.0001)
# model.compile(loss=['categorical_crossentr'], optimizer=opt)

# # Training


# history1 = model.fit(x_train, y_train, epochs=epochsToUse, batch_size=batchSizeToUse, validation_data=(x_validate, y_validate))

# #Plotting
# plt.plot(history1.history['mse'])
# plt.plot(history1.history['val_mse'])
# plt.title('Accuracy vs Epochs for Classifier')
# plt.ylabel('accuracy')
# plt.xlabel('epoch')
# plt.legend(['train', 'test'], loc='upper left')
# plt.show()

plt.plot(history1.history['loss'])
plt.plot(history1.history['val_loss'])
plt.title('Loss vs Epochs for Classifier')
plt.ylabel('loss')
plt.xlabel('epoch')
plt.legend(['train', 'test'], loc='upper left')
plt.show()


# plt.plot(history1.history['loss'])
# plt.plot(history1.history['val_loss'])
# plt.title('Loss vs Epochs for Classifier')
# plt.ylabel('loss')
# plt.xlabel('epoch')
# plt.legend(['train', 'test'], loc='upper left')
# plt.show()







vectorImagesTest = np.zeros((0, sizeOfImg, sizeOfImg, 3))
outputLayerTest = np.zeros((0, 7))


for base, dirs, files in os.walk(test_folder_path):
    # Iterate over the files
    #print("We are here: ",type(files))
    files.sort()
    
    for file in files:
        # Check if the file has a JPEG or JPG extension
        #print(file,'\n')
        
        if word_to_find.lower() in file.lower():
            print("not using files: ",word_to_find)
        else:
            if file.lower().endswith(('.jpeg', '.jpg')):
                # If yes, construct the full file path and add it to the list
                file_path = os.path.join(base, file)
                dslocal, original_size = resize_image(file_path, sizeOfImg, sizeOfImg)
                # print(dslocal.shape)
                vectorImagesTest = np.append(vectorImagesTest, [dslocal], axis=0)
                # print(list(dslocal.shape))
                currentX = original_size[0]
                currentY = original_size[1]
                
                print(file_path)
            elif file.lower().endswith(('.xml')):
                labels, bboxes = parse_xml(test_folder_path+file)
                # print(labels)
                # print(bboxes)
                
                
                
                bndWidth = int( ( ( list(bboxes[0])[2] - list(bboxes[0])[0] )/currentX )*sizeOfImg )
                bndHeight = int( ( ( list(bboxes[0])[3] - list(bboxes[0])[1] )/currentY )*sizeOfImg )
                xCenter = int( ( ( list(bboxes[0])[2] + list(bboxes[0])[0] )/ (currentX*2) )*sizeOfImg )
                yCenter = int( ( ( list(bboxes[0])[3] + list(bboxes[0])[1] )/ (currentY*2) )*sizeOfImg )
                
                # bndWidth = int( ( ( list(bboxes[0])[0] )/currentX )*sizeOfImg )
                # bndHeight = int( ( ( list(bboxes[0])[1] )/currentY )*sizeOfImg )
                # xCenter = int( ( ( list(bboxes[0])[2] )/ (currentX) )*sizeOfImg )
                # yCenter = int( ( ( list(bboxes[0])[3] )/ (currentY) )*sizeOfImg )
                
                
                
                if labels[0] == "apple":
                    localArray = np.array([1, 0, 0, xCenter, yCenter, bndWidth, bndHeight])
                    print(localArray)
                    outputLayerTest = np.append(outputLayerTest, [localArray], axis=0)
                elif labels[0] == "banana":
                    localArray = np.array([0, 1, 0, xCenter, yCenter, bndWidth, bndHeight])
                    print(localArray)
                    outputLayerTest = np.append(outputLayerTest, [localArray], axis=0)
                    # outputLayerTest.append([1, 1, list(bboxes[0])])
                elif labels[0] == "orange":
                    localArray = np.array([0, 0, 1, xCenter, yCenter, bndWidth, bndHeight])
                    print(localArray)
                    outputLayerTest = np.append(outputLayerTest, [localArray], axis=0)
                    # outputLayerTest.append([1, 2, list(bboxes[0])])
       
            

print("Testing X shape: "+str(vectorImagesTest.shape))
print("Testing Y shape: "+str(outputLayerTest.shape))



predictions = model.predict(vectorImagesTest)
# predicted_labels = np.argmax(predictions, axis=1)
# true_labels = outputLayerTest
# accuracy = accuracy_score(true_labels, predicted_labels)
# report = classification_report(true_labels, predicted_labels)

# print("Accuracy:", accuracy)
# print("Classification Report:\n", report)












