#!/usr/bin/env python
# coding: utf-8

# # Import Dataset

# In[1]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter
from sklearn.preprocessing import LabelBinarizer
from keras.preprocessing import text, sequence
from keras.preprocessing.text import Tokenizer
from sklearn.model_selection import train_test_split
from keras import backend as K
import tensorflow as tf
from keras.models import Sequential
from keras.layers import Dense, Conv1D, MaxPooling1D, Flatten
from keras.layers import LSTM, GRU, Bidirectional, Dropout, SimpleRNN, Input
from keras.layers import Embedding
from keras.models import Model
from keras.callbacks import TensorBoard
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import precision_recall_fscore_support
from keras.utils import to_categorical
from keras.callbacks import EarlyStopping
from sklearn.metrics import confusion_matrix
from sklearn.metrics import accuracy_score
from sklearn.metrics import classification_report
import itertools
from keras.layers import BatchNormalization
import os
from keras.engine.topology import Layer

def dataset(top_classes = 4, path = 'full_pair_data.csv'):
    df = pd.read_csv(path)
    df.reset_index()
    cnt = Counter(df.Class)
    sorted_classes = cnt.most_common()[:int(top_classes)]
    classes = [c[0] for c in sorted_classes]
    counts = [c[1] for c in sorted_classes] 
    df = df[[c in classes for c in df.Class]]
    print(str(df.shape[0]) + " instances")
    seqs = df.Seqs.values
    lengths = [len(s) for s in seqs]
    lb = LabelBinarizer()
    Y = lb.fit_transform(df.Class)
    max_length = 512
    seqs_sys = 'ARNDCQEGHILKMFPSTWYV-'
    tokenizer = Tokenizer(char_level=True)
    tokenizer.fit_on_texts(seqs_sys)
    X = tokenizer.texts_to_sequences(seqs)
    X = sequence.pad_sequences(X, maxlen = max_length)   
    return X, Y


class SelfAttention(Layer):
    def __init__(self, output_dim, **kwargs):
        self.output_dim = output_dim
        super(SelfAttention, self).__init__(**kwargs)

    def build(self, input_shape):
        self.kernel = self.add_weight(name='kernel',
                                      shape=(3, input_shape[2], self.output_dim),
                                      initializer = 'uniform',
                                      trainable = True)
        super(SelfAttention, self).build(input_shape)

    def get_config(self):
        config = super(SelfAttention, self).get_config().copy()
        config.update({
            'output_dim': self.output_dim,
        })
        return config

    def call(self, x):
        WQ = K.dot(x, self.kernel[0])
        WK = K.dot(x, self.kernel[1])
        WV = K.dot(x, self.kernel[2])
        # print("WQ.shape",WQ.shape)
        # print("K.permute_dimensions(WK, [0, 2, 1]).shape",K.permute_dimensions(WK, [0, 2, 1]).shape)
        QK = K.batch_dot(WQ,K.permute_dimensions(WK, [0, 2, 1]))
        QK = QK / (self.output_dim**0.5)
        QK = K.softmax(QK)
        # print("QK.shape",QK.shape)
        return K.batch_dot(QK,WV)

    def compute_output_shape(self, input_shape):
        return (input_shape[0],input_shape[1],self.output_dim)



class DLPAlign():
    def __init__(self, mode = 'train', model_addr = 'dlpalign.h5', top_classes = 4, input_length = 512, embedding_dim = 8, data_path = 'result.csv'):
        self.input_length = input_length
        self.embedding_dim = embedding_dim
        self.data_path = data_path
        self.top_classes = top_classes
        self.input_dim = 22
        self.X, self.Y = dataset(self.top_classes, self.data_path)
        self.model = self.structure(self.input_length, self.embedding_dim)
        if mode == 'test':
            self.model.load_weights(model_addr)
        self.model.compile(loss = self.multi_category_focal_loss, optimizer = 'adam', metrics=['accuracy'])

    def multi_category_focal_loss(self, y_true, y_pred):
        epsilon = 1.e-7
        gamma = 2.0
        alpha = tf.constant([[2],[3],[2],[1]], dtype=tf.float32)
        y_true = tf.cast(y_true, tf.float32)
        y_pred = tf.clip_by_value(y_pred, epsilon, 1. - epsilon)
        y_t = tf.multiply(y_true, y_pred) + tf.multiply(1-y_true, 1-y_pred)
        ce = -tf.math.log(y_t)
        weight = tf.pow(tf.subtract(1., y_t), gamma)
        fl = tf.matmul(tf.multiply(weight, ce), alpha)
        loss = tf.reduce_mean(fl)
        return loss

    def structure(self, input_length, embedding_dim):
        inputs = Input(shape=(input_length,))
        l1 = Embedding(int(self.input_dim), int(embedding_dim), input_length = int(input_length))(inputs)
        l2 = Conv1D(filters=32, kernel_size=3, padding='same', activation='selu')(l1)
        l2 = BatchNormalization()(l2)
        l2 = MaxPooling1D(pool_size=2)(l2)
        l3 = Conv1D(filters=64, kernel_size=3, padding='same', activation='selu')(l2)
        l3 = BatchNormalization()(l3)
        l3 = MaxPooling1D(pool_size=2)(l3)
        l4 = Conv1D(filters=128, kernel_size=3, padding='same', activation='selu')(l3)
        l4 = BatchNormalization()(l4)
        l4 = MaxPooling1D(pool_size=2)(l4)
        l5 = Conv1D(filters=32, kernel_size=3, padding='same', activation='selu')(l4)
        l5 = BatchNormalization()(l5)
        l5 = MaxPooling1D(pool_size=2)(l5)
        l6 = Bidirectional(LSTM(128, return_sequences=True))(l5)
        l6 = BatchNormalization()(l6)
        atten = SelfAttention(64)(l6)
        atten = Flatten()(atten)
        l7 = Dense(256, activation='selu')(atten)
        l7 = Dropout(rate=0.5)(l7)
        l8 = Dense(32, activation='selu')(l7)
        pred = Dense(int(self.top_classes), activation='softmax')(l8)
        model = Model(inputs = inputs, outputs = pred)
        model.summary()
        return model

    def train(self, epochs = 30, batch_size = 128, save_path = 'dlpalign.h5'):
        early_stopping = EarlyStopping(monitor='val_loss', patience=20, mode='auto')
        seed = 929
        np.random.seed(seed)
        kfold = StratifiedKFold(n_splits = 4, shuffle=True, random_state=seed)
        callbacks = [early_stopping]
        X_train, X_test, y_train, y_test = train_test_split(self.X, self.Y.argmax(1), test_size=.2)
        for train, test in kfold.split(X_train, y_train):
            self.model.fit(X_train[train], to_categorical(y_train[train]), validation_data = (X_train[test], to_categorical(y_train[test])), epochs = epochs, batch_size = batch_size, callbacks = callbacks)
        self.test(X_test, y_test)
        self.model.save(save_path)

    def test(self, X_test, y_test):
        print(len(X_test))
        test_pred = self.model.predict(X_test)
        print("test-acc = " + str(accuracy_score(y_test, np.argmax(test_pred, axis=1))))
        y_pred = np.argmax(test_pred, axis=1)
        with open('y_test_pred.txt', 'w') as fileout:
            for i in range(len(y_test)):
                fileout.write(str(i) + '\t' + str(y_test[i]) + '\t' + str(y_pred[i]) + '\n')
        print('Done!')
        print(classification_report(y_test.tolist(), np.argmax(test_pred, axis = 1).tolist(), digits = 4, target_names=['0', '1', '2', '3']))

if __name__ == '__main__':
    dlpalign = DLPAlign()
    dlpalign.train()