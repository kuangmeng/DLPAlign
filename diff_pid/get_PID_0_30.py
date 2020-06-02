#!/usr/bin/env python
import  xdrlib ,sys
import xlrd
from xlwt import *
import math, os

txt_file = ""

def open_excel(file):
    data = xlrd.open_workbook(file)
    return data

def ReadTXT(txt_file):
    ret_matrix30 = []
    with open(txt_file, 'r') as file_in:
        file_context = file_in.read().splitlines()
        for line in file_context:
            tmp_list = line.split(", ")
            if float(tmp_list[1]) <= 0.3:
                ret_matrix30.append(tmp_list[0])
    return ret_matrix30

def ReadEXCEL(excel_file, ret_matrix):
    ret_1 = 0.0
    ret_2 = 0.0
    lens = len(ret_matrix)
    data = open_excel(excel_file)
    table = data.sheets()[0]
    nrows = table.nrows
    best_ = 0.0
    best = [0 for i in range(2)]
    for line in ret_matrix:
        for j in range(nrows):
            row = table.row_values(j)
            if line == row[0]:
                ret_1 += float(row[1])
                ret_2 += float(row[2])
                break
    return ret_1 / lens, ret_2 / lens


if __name__ == '__main__':
    txt_file = "pid.txt"
    file_list = os.listdir('./')
    ret_1 = ReadTXT(txt_file)
    print(len(ret_1))
    print("PID < 30:")
    for file in file_list:
        if '.xlsx' in file:
            excel = file
            tc, sp = ReadEXCEL(excel, ret_1)
            print("%s: %f" %(excel, tc))
