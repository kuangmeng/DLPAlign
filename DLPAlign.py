#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 31 16:06:43 2020

@author: kuangmeng
"""
import os
import subprocess
import sys
from joblib import load
from pair_testing import SingFileTest
from tensorflow.keras.models import load_model
import time

first_step = 0
second_step = 0
third_step = 0

model = load_model('./Classifier/model/pairs_classification_cnn_bilstm_final4.h5')

def ExistOrNot():
    return not os.path.exists('dlpalign_main')

def Clean_Make():
    log = open('./log.tmp', 'w')
    _, clean = subprocess.getstatusoutput('make clean')
    log.write('<-- make clean log-->\n')
    log.write(clean + '\n')
    _, make = subprocess.getstatusoutput('make')
    log.write('<-- make log -->\n')
    log.write(make + '\n')
    log.close()

# def Normalization(tmp_list, para):
#     paras = []
#     with open(para, 'r') as filein:
#         paras = filein.read().splitlines()
#     for i in range(len(tmp_list)):
#         tmp_list[i] = (float(tmp_list[i]) - float(paras[2 * i + 1])) / (float(paras[i * 2]) - float(paras[i * 2 - 1]))
#     return [tmp_list]

# def getTestList(in_file):
#     _, tmp_str = subprocess.getstatusoutput('./mlpalign_main -G ' + in_file)
#     test_list = tmp_str.split(',')
#     return test_list

def getClassification(in_file, model):
    result = SingFileTest(in_file, model)
    return int(result)

def MakeCMD(first_step, second_step, third_step, in_file, out_file):
    cmd = './dlpalign_main -s1 ' + str(first_step) + ' -s2 ' + str(second_step) + ' -s3 ' + str(third_step) + ' ' + in_file
    os.system(cmd + ' > ' + out_file)
    print('DLPAlign Finished!')

out_file = 'result.msa'

def RunSingle(in_file, out_file, model):
    print('Input: %s'%(in_file))
    start = time.time()
    first_step = getClassification(in_file, model)
    MakeCMD(first_step, second_step, third_step, in_file, out_file)
    end = time.time()
    return float(end) - float(start)
    
def RunMultiple(in_bench, out_dir, model):
    time = 0
    for file in os.listdir(in_bench + '/in' ):
        if file[0] != '.':
            in_file = in_bench + '/in/' + file
            out_file = out_dir + '/' + in_bench.split('/')[-1] + '/' + file
            if not os.path.exists(out_dir + '/' + in_bench.split('/')[-1]):
                os.system('mkdir %s'%(out_dir + '/' + in_bench.split('/')[-1]))
            time += RunSingle(in_file, out_file, model)
    print('%s: %f'%(in_bench, time))

if __name__ == '__main__':
    # if ExistOrNot():
    #     Clean_Make()
    out_dir = './output'
    if not os.path.exists(out_dir):
                os.system('mkdir %s'%(out_dir))
    in_bench = sys.argv[1]
    RunMultiple(in_bench, out_dir, model)
    
    
    
        
    
        



