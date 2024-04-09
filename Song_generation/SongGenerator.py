#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


@author: josue_nataren
"""

#Song generator

import tensorflow as tf
from tensorflow import keras
from tensorflow.keras.preprocessing.text import Tokenizer
from tensorflow.keras.preprocessing.sequence import pad_sequences
import numpy as np


tokenizer = Tokenizer()

with open('AdeleBadBunny21Pilots.txt', 'r') as file:
    # Read the entire contents of the file into a variable
    fileContent = file.read()

# Now file_content variable contains all the content of the file
print(fileContent)

corpusToUse = fileContent.lower().split("\n")

tokenizer.fit_on_texts(corpusToUse)

totalNumber = len(tokenizer.word_index) + 1


print(totalNumber)


inputSequence = []

for line in corpusToUse:
    tokenList = tokenizer.texts_to_sequences([line])[0]
    for i in range(1, len(tokenList)):
        nGram_sequence = tokenList[:i+1]
        inputSequence.append(nGram_sequence)


maxSequenceLength = max([len(x) for x in inputSequence])

inputSequence = np.array(pad_sequences(inputSequence, maxlen = maxSequenceLength, padding="pre"))


xs = inputSequence[:,:-1]
labels = inputSequence[:,-1]

ys = tf.keras.utils.to_categorical(labels, num_classes = totalNumber)


model = tf.keras.Sequential()
model.add(tf.keras.layers.Embedding(totalNumber, 180, input_length=maxSequenceLength-1))
model.add(tf.keras.layers.Bidirectional(tf.keras.layers.LSTM(256, return_sequences=True)))
model.add(tf.keras.layers.Dropout(0.2))
model.add(tf.keras.layers.Bidirectional(tf.keras.layers.LSTM(128)))
model.add(tf.keras.layers.Dropout(0.2))
model.add(tf.keras.layers.Dense(totalNumber, activation = 'softmax'))

adam = tf.keras.optimizers.Adam(learning_rate = .001)

model.compile( loss = 'categorical_crossentropy', optimizer = adam, metrics=['accuracy'] )

history = model.fit(xs, ys, epochs = 70, verbose = 1)


initialText = "Nice to know"
nextWords = 120


for _ in range(nextWords):
    tokenList = tokenizer.texts_to_sequences([initialText])[0]
    tokenList = pad_sequences([tokenList], maxlen = maxSequenceLength-1, padding="pre")
    predicted = model.predict(tokenList) 
    classesX = np.argmax(predicted,axis=1)
    nextLocalWord = ""
    for word, index in tokenizer.word_index.items():
        if index == classesX:
            nextLocalWord = word
            break
        
    initialText += " " + nextLocalWord
    

print(initialText)





