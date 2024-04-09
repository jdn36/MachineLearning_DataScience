#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: josue_nataren
"""

import numpy as np
import matplotlib.pyplot as plt
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense
import random
from sklearn.metrics import accuracy_score, classification_report
from sklearn.model_selection import train_test_split
import pywt
import tensorflow as tf
from tensorflow.keras.utils import to_categorical

train_raw = np.genfromtxt('MIT-BIH/mitbih_train.csv',delimiter=',')
test_raw = np.genfromtxt('MIT-BIH/mitbih_test.csv',delimiter=',')


[rows, cols] = train_raw.shape
[rowsTest, colsTest] = test_raw.shape

#Check plot
plt.plot(train_raw[0][:cols-1])
plt.show()

x_train = train_raw[:rows, :cols-1]
y_train = train_raw[:rows, cols-1]

x_test = test_raw[:rowsTest, :colsTest-1]
y_test = test_raw[:rowsTest, colsTest-1]

level_Decomposition = 2
x_coeffs_train = []
x_coeffs_test = []

for i in range(rows):
    #x_coeffs_train.append(pywt.wavedec(x_train[i], 'sym5', level=level_Decomposition))
    cA2, cD2, cD1 = pywt.wavedec(x_train[i], 'sym5', level=level_Decomposition)
    x_coeffs_train.append(cA2)

for i in range(rowsTest):
    #x_coeffs_train.append(pywt.wavedec(x_train[i], 'sym5', level=level_Decomposition))
    cA2, cD2, cD1 = pywt.wavedec(x_test[i], 'sym5', level=level_Decomposition)
    x_coeffs_test.append(cA2)


x_coeffs_train = np.array(x_coeffs_train)
x_coeffs_test = np.array(x_coeffs_test)

print(x_coeffs_train.shape)
print(x_coeffs_test.shape)


plt.plot(x_coeffs_train[0, :53])
plt.show


num_Classes = np.unique(y_train).size
[rowsInput, colsInput] = x_coeffs_train.shape

y_train_matrix = to_categorical(y_train - 1, num_classes=num_Classes)
y_test_matrix = to_categorical(y_test - 1, num_classes=num_Classes)

x_train_coeffs, x_test_coeffs, y_train_coeffs, y_test_coeffs = train_test_split(x_coeffs_train, y_train_matrix, test_size=0.2, random_state=42, stratify=y_train_matrix)

# Model architecture
model = Sequential([
    Dense(128, activation='relu', input_shape=(colsInput,)),
    Dense(256, activation='relu'),
    Dense(512, activation='relu'),
    Dense(256, activation='relu'),
    Dense(128, activation='relu'),
    Dense(64, activation='relu'),
    Dense(32, activation='relu'),
    Dense(num_Classes, activation='softmax')
])

model.compile(optimizer='adam',
              loss='categorical_crossentropy',
              metrics=['accuracy'])

model.summary()

epochs_Num = 25
batch_size_Num = 20
history1 = model.fit(x_train_coeffs, y_train_coeffs, epochs=epochs_Num, batch_size=batch_size_Num, validation_data=(x_coeffs_test, y_test_matrix))


#Plotting
plt.plot(history1.history['accuracy'])
plt.plot(history1.history['val_accuracy'])
plt.title('Accuracy vs Epochs for Classifier')
plt.ylabel('accuracy')
plt.xlabel('epoch')
plt.legend(['train', 'test'], loc='upper left')
plt.show()

plt.plot(history1.history['loss'])
plt.plot(history1.history['val_loss'])
plt.title('Loss vs Epochs for Classifier')
plt.ylabel('loss')
plt.xlabel('epoch')
plt.legend(['train', 'test'], loc='upper left')
plt.show()


predictions = model.predict(x_test_coeffs)
predicted_labels = np.argmax(predictions, axis=1)
true_labels = np.argmax(y_test_coeffs, axis=1)
accuracy = accuracy_score(true_labels, predicted_labels)
report = classification_report(true_labels, predicted_labels)

print("Accuracy:", accuracy)
print("Classification Report:\n", report)


if history1.history['accuracy'][-1] > .9:
    print('Training accuracy higher than 90%')

if history1.history['val_accuracy'][-1] > .9:
    print('Validation accuracy higher than 90%')

if accuracy > .9:
    print('Testing accuracy higher than 90%')


# array1 = np.array(history1.history['accuracy'])
# array2 = np.array(history1.history['val_accuracy'])
# column_vector1 = array1.reshape(-1, 1)
# column_vector2 = array2.reshape(-1, 1)
# exportingAccuracy = np.concatenate((column_vector1, column_vector2), axis=1)
# np.savetxt('modelAccuracy2.csv', exportingAccuracy, delimiter=',')

# model.save_weights('model_weights2.h5')
# model.save('full_model2.h5')
# model.save('full_model2.keras')





