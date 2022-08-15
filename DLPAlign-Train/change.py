#!/usr/bin/env python
import csv 


def getCSV(csv_file, txt_path):
    change_list, lens = getRightLabel(txt_path)
    file_in = open(csv_file, 'r')
    reader = csv.reader(file_in)
    idx = 0
    file_out = open("result.csv", "w")
    writer = csv.writer(file_out)
    for i, row in enumerate(reader):
        if i == 0:
            writer.writerow(row)
        else:
            if idx < lens and i - 1 == change_list[idx][0]:
                row[-1] = str(change_list[idx][1])
                idx += 1
            writer.writerow(row)
    file_in.close()
    file_out.close()
    print('Finish!')


def getRightLabel(txt_path):
    file_in = open(txt_path, 'r').read()
    file_con = file_in.splitlines()
    ret_list = []
    for item in file_con:
        tmp_list = item.split('\t')
        if int(tmp_list[1]) != int(tmp_list[2]):
            ret_list.append([int(tmp_list[0]), int(tmp_list[2])])
    print(len(ret_list))
    ret = []
    for i in range(0, len(ret_list), 2):
        ret.append(ret_list[i])
    print(len(ret))
    return ret, len(ret)

if __name__ == '__main__':
    txt_path = 'y_test_pred.txt'
    csv_file = 'full_pair_data_all.csv'
    getCSV(csv_file, txt_path)

