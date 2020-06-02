#!/usr/bin/env python
# coding: utf-8

# # Import Dataset

# In[1]:


import numpy as np
import pandas as pd

df = pd.read_csv('./input/full_pair_data_3.csv')
df.reset_index()
df.shape


# # Preprocessing and visualization of dataset

# In[2]:


import matplotlib.pyplot as plt
from collections import Counter

# count numbers of instances per class
cnt = Counter(df.Class)
# select only 10 most common classes!
top_classes = 4
# sort classes
sorted_classes = cnt.most_common()[:top_classes]
classes = [c[0] for c in sorted_classes]
counts = [c[1] for c in sorted_classes]
print(classes)
print(counts)

print("at least " + str(counts[-1]) + " instances per class")

# apply to dataframe
print(str(df.shape[0]) + " instances before")
df = df[[c in classes for c in df.Class]]
print(str(df.shape[0]) + " instances after")

seqs = df.Seqs.values
lengths = [len(s) for s in seqs]

# # visualize
# fig, axarr = plt.subplots(1,2, figsize=(20,5))
# axarr[0].bar(range(len(classes)), counts)
# plt.sca(axarr[0])
# plt.xticks(range(len(classes)), classes, rotation='vertical')
# axarr[0].set_ylabel('frequency')

# axarr[1].hist(lengths, bins=100, normed=False)
# axarr[1].set_xlabel('sequence length')
# axarr[1].set_ylabel('# sequences')
# plt.show()


# # Transform labels

# In[3]:


from sklearn.preprocessing import LabelBinarizer

# Transform labels to one-hot
lb = LabelBinarizer()
Y = lb.fit_transform(df.Class)


# # Preprocessing 

# In[4]:


from keras.preprocessing import text, sequence
from keras.preprocessing.text import Tokenizer
from sklearn.model_selection import train_test_split

# maximum length of sequence, everything afterwards is discarded!
max_length = 512

seqs_sys = 'ARNDCQEGHILKMFPSTWYV-'
#create and fit tokenizer
tokenizer = Tokenizer(char_level=True)
tokenizer.fit_on_texts(seqs_sys)
#represent input data as word rank number sequences
X = tokenizer.texts_to_sequences(seqs)
X = sequence.pad_sequences(X, maxlen=max_length)


# # Build keras model and fit

# In[5]:


from keras import backend as K
import tensorflow as tf

def multi_category_focal_loss(y_true, y_pred):
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


# In[14]:


from  tensorflow.keras.models import Sequential
from  tensorflow.keras.layers import Dense, Conv1D, MaxPooling1D, Flatten
from  tensorflow.keras.layers import LSTM, GRU, Bidirectional, Dropout, SimpleRNN
from  tensorflow.keras.layers import Embedding
from  keras.utils.vis_utils import plot_model

embedding_dim = 8

# model = Sequential()
# model.add(Embedding(len(tokenizer.word_index)+1, embedding_dim, input_length=max_length))
# model.add(Conv1D(filters=64, kernel_size=6, padding='same', activation='selu'))
# model.add(MaxPooling1D(pool_size=2))
# model.add(Conv1D(filters=32, kernel_size=3, padding='same', activation='selu'))
# model.add(MaxPooling1D(pool_size=2))
# model.add(Bidirectional(LSTM(64, return_sequences=True)))
# #model.add(GRU(64, activation='tanh', recurrent_activation='hard_sigmoid', use_bias=True, kernel_initializer='glorot_uniform', recurrent_initializer='orthogonal', bias_initializer='zeros', kernel_regularizer=None, recurrent_regularizer=None, bias_regularizer=None, activity_regularizer=None, kernel_constraint=None, recurrent_constraint=None, bias_constraint=None, dropout=0.0, recurrent_dropout=0.0, implementation=1, return_sequences=True, return_state=False, go_backwards=False, stateful=False, unroll=False, reset_after=False))
# #model.add(SimpleRNN(16, return_sequences=True, activation='selu'))
# #model.add(LSTM(64, return_sequences=True))
# model.add(Flatten())
# model.add(Dense(1024, activation='selu'))
# model.add(Dropout(rate=0.5))
# model.add(Dense(128, activation='selu'))
# model.add(Dense(top_classes, activation='softmax'))
# model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
# print(model.summary())

model = Sequential()
model.add(Embedding(len(tokenizer.word_index)+1, embedding_dim, input_length=max_length))
model.add(Conv1D(filters=32, kernel_size=6, padding='same', activation='selu'))
model.add(MaxPooling1D(pool_size=2))
model.add(Conv1D(filters=16, kernel_size=3, padding='same', activation='selu'))
model.add(MaxPooling1D(pool_size=2))
#model.add(Bidirectional(GRU(16, return_sequences=True)))
#model.add(GRU(64, activation='tanh', recurrent_activation='hard_sigmoid', use_bias=True, kernel_initializer='glorot_uniform', recurrent_initializer='orthogonal', bias_initializer='zeros', kernel_regularizer=None, recurrent_regularizer=None, bias_regularizer=None, activity_regularizer=None, kernel_constraint=None, recurrent_constraint=None, bias_constraint=None, dropout=0.0, recurrent_dropout=0.0, implementation=1, return_sequences=True, return_state=False, go_backwards=False, stateful=False, unroll=False, reset_after=False))
model.add(SimpleRNN(16, return_sequences=True, activation='selu'))
#model.add(GRU(16, return_sequences=True))
model.add(Flatten())
model.add(Dense(256, activation='selu'))
model.add(Dropout(rate=0.5))
model.add(Dense(64, activation='selu'))
model.add(Dense(top_classes, activation='softmax'))
model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
print(model.summary())


# In[15]:


from tensorflow.keras.callbacks import TensorBoard
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import precision_recall_fscore_support
from keras.utils import to_categorical
from tensorflow.keras.callbacks import EarlyStopping

early_stopping = EarlyStopping(monitor='val_loss', patience=10, mode='auto')

seed = 929
np.random.seed(seed)
kfold = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
callbacks = [early_stopping]

X_train, X_test, y_train, y_test = train_test_split(X, Y.argmax(1), test_size=.2)

for train, test in kfold.split(X_train, y_train):
    model.fit(X_train[train], to_categorical(y_train[train]), validation_data=(X_train[test], to_categorical(y_train[test])), epochs=20, batch_size=128, callbacks=callbacks)

#model.save("pairs_classification_cnn_bilstm.h5")



# # Evaluate

# In[16]:


from sklearn.metrics import confusion_matrix
from sklearn.metrics import accuracy_score
from sklearn.metrics import classification_report
import itertools

#train_pred = model.predict(X_train)
test_pred = model.predict(X_test)
#print("train-acc = " + str(accuracy_score(np.argmax(y_train, axis=1), np.argmax(train_pred, axis=1))))
print("test-acc = " + str(accuracy_score(y_test, np.argmax(test_pred, axis=1))))

y_pred = np.argmax(test_pred, axis=1)

with open('y_test_pred.txt', 'w') as fileout:
    for i in range(len(y_test)):
        fileout.write(str(y_test[i]) + '\t' + str(y_pred[i]) + '\n')
print('Done!')
# # Compute confusion matrix
# cm = confusion_matrix(np.argmax(y_test, axis=1), np.argmax(test_pred, axis=1))

# # Plot normalized confusion matrix
# cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
# np.set_printoptions(precision=2)
# plt.figure(figsize=(10,10))
# plt.imshow(cm, interpolation='nearest', cmap=plt.cm.Blues)
# plt.title('Confusion matrix')
# plt.colorbar()
# tick_marks = np.arange(len(lb.classes_))
# plt.xticks(tick_marks, lb.classes_, rotation=90)
# plt.yticks(tick_marks, lb.classes_)
# #for i, j in itertools.product(range(cm.shape[0]), range(cm.shape[1])):
# #    plt.text(j, i, format(cm[i, j], '.2f'), horizontalalignment="center", color="white" if cm[i, j] > cm.max() / 2. else "black")
# plt.ylabel('True label')
# plt.xlabel('Predicted label')
# plt.show()


# In[17]:


print(classification_report(y_test.tolist(), np.argmax(test_pred, axis=1).tolist(), digits=4, target_names=['0', '1', '2', '3']))


# In[ ]:




