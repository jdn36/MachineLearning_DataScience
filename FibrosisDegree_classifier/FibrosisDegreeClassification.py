import matplotlib.pyplot as plt
import numpy as np
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense
import random
from sklearn.metrics import accuracy_score, classification_report
# import csv

matrixTraining = np.genfromtxt('./fibrosisDegreeTraining.csv', delimiter=',')
matrixTest = np.genfromtxt('./fibrosisDegreeTest.csv', delimiter=',')
matrixOther= np.genfromtxt('./fibrosisDegreeOther.csv', delimiter=',')

vectorTraining = np.genfromtxt('./fibrosisDegreeTrainingAns.csv', delimiter=',')
vectorTest = np.genfromtxt('./fibrosisDegreeTestAns.csv', delimiter=',')
vectorOther= np.genfromtxt('./fibrosisDegreeOtherAns.csv', delimiter=',')
        
print(type(matrixTraining))
print(matrixTraining.shape)
# print(matrixTraining[num_rows_Training-2][num_columns_Training-2])

num_rows_Training, num_columns_Training = matrixTraining.shape
num_rows_Test, num_columns_Test = matrixTest.shape

num_rows_Outpus, num_columns_Ouput = vectorTraining.shape

indicesToUseTraining = [i for i in range(num_rows_Training)]
random.shuffle(indicesToUseTraining)
# print(indicesToUseTraining)
matrixTraining = matrixTraining[indicesToUseTraining][:]
vectorTraining = vectorTraining[indicesToUseTraining][:]

indicesToUseTest = [i for i in range(num_rows_Test)]
random.shuffle(indicesToUseTest)
# print(indicesToUseTest)
matrixTest = matrixTest[indicesToUseTest][:]
vectorTest = vectorTest[indicesToUseTest][:]


# Model architecture
model = Sequential([
    Dense(16, activation='relu', input_shape=(num_columns_Training,)),
    Dense(32, activation='relu'),
    Dense(64, activation='relu'),
    Dense(128, activation='relu'),
    Dense(64, activation='relu'),
    Dense(32, activation='relu'),
    Dense(num_columns_Ouput, activation='softmax')
])



model.compile(optimizer='adam',
              loss='categorical_crossentropy',
              metrics=['accuracy'])

model.summary()

epochs_Num = 90
batch_size_Num = 16
history1 = model.fit(matrixTraining, vectorTraining, epochs=epochs_Num, batch_size=batch_size_Num, validation_data=(matrixTest, vectorTest))


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


predictions = model.predict(matrixOther)
predicted_labels = np.argmax(predictions, axis=1)
true_labels = np.argmax(vectorOther, axis=1)
accuracy = accuracy_score(true_labels, predicted_labels)
report = classification_report(true_labels, predicted_labels)

print("Accuracy:", accuracy)
print("Classification Report:\n", report)





array1 = np.array(history1.history['accuracy'])
array2 = np.array(history1.history['val_accuracy'])
column_vector1 = array1.reshape(-1, 1)
column_vector2 = array2.reshape(-1, 1)
exportingAccuracy = np.concatenate((column_vector1, column_vector2), axis=1)
np.savetxt('modelAccuracy2.csv', exportingAccuracy, delimiter=',')

array1 = np.array(history1.history['loss'])
array2 = np.array(history1.history['val_loss'])
column_vector1 = array1.reshape(-1, 1)
column_vector2 = array2.reshape(-1, 1)
exportingAccuracy = np.concatenate((column_vector1, column_vector2), axis=1)
np.savetxt('modelLoss2.csv', exportingAccuracy, delimiter=',')


model.save_weights('model_weights2.h5')
model.save('full_model2.h5')
model.save('full_model2.keras')



