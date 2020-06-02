from tensorflow.keras.preprocessing import text, sequence
from tensorflow.keras.preprocessing.text import Tokenizer
from sklearn.model_selection import train_test_split
from tensorflow.keras.models import load_model
import numpy as np
import os, sys
from pandas.core.frame import DataFrame

max_length = 512
seqs_sys = 'ARNDCQEGHILKMFPSTWYV-'

def FindMax(np_arr, lens):
    list_sys = [0.22, 0.17, 0.21, 0.4]
    max = 0
    max_idx = 0
    for i in range(len(np_arr)):
        if (np_arr[i] / lens) / list_sys[i] > max:
            max_idx = i
            max = (np_arr[i] / lens) / list_sys[i]
    return max_idx


def SingleTest(seqs, model):
    df = DataFrame({'Seqs': seqs})
    df.reset_index()
    seqs_ = df.Seqs.values
    tokenizer = Tokenizer(char_level=True)
    tokenizer.fit_on_texts(seqs_sys)
    X = tokenizer.texts_to_sequences(seqs_)
    X = sequence.pad_sequences(X, maxlen=max_length)
    np_list = np.argmax(model.predict(X), axis=1)
    return FindMax(np.bincount(np_list), len(np_list))

def SingFileTest(family_file, model):
    seqs = []
    tmp_val = ""
    pairs = []
    filein = open(family_file, 'r')
    file_con = filein.read().splitlines()
    filein.close()
    for line in file_con:
        if len(tmp_val) > 0 and '>' in line:
            seqs.append(tmp_val)
            tmp_val = ""
        elif '>' not in line:
            tmp_val += line
    if len(tmp_val) > 0:
        seqs.append(tmp_val)
    for i in range(len(seqs) - 1):
        for j in range(i, len(seqs)):
            pairs.append(seqs[i] + '-' + seqs[j])
    return SingleTest(pairs, model)

def BenchTest(bench, model):
    ret = []
    files = os.listdir(bench + '/in')
    for file in files:
        if file[0] != '.':
            ret.append([bench + '.' + file, SingFileTest(bench + '/in/' +  file, model)])
    return ret
    

if  __name__ == "__main__":
    bench = sys.argv[1]
    ret = BenchTest(bench, model)
    print(ret)