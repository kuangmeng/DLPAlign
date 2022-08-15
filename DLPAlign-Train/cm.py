#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 21 19:02:52 2021

@author: kuangmeng
"""
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay
import numpy as np


class CM:
    def __init__(self, txt_path, classes = [0, 1, 2, 3]):
        self.y_true, self.pred = self.parser(txt_path)
        self.classes = classes
        
    def parser(self, txt_path):
        y_true = []
        pred = []
        with open(txt_path, 'r') as filein:
            file_con = filein.readlines()
            for item in file_con:
                y_true.append(int(item.split()[1]))
                pred.append(int(item.split()[2]))
        return np.array(y_true), np.array(pred)
    
    def dispaly(self):
        cm = confusion_matrix(self.y_true, self.pred, labels = self.classes)
        disp = ConfusionMatrixDisplay(confusion_matrix = cm, display_labels = self.classes)
        disp.plot() 
        
        
        
if __name__ == '__main__':
    cm = CM('y_test_pred.txt')
    cm.dispaly()