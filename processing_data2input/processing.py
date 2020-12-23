import csv
import os.path
from os import path


# create file from the data that the user inserted manually to web
def add2file(file_path, famcode, A1_data, B1_data, C1_data, DRB11_data, DQB11_data, A2_data, B2_data, C2_data, DRB12_data, DQB12_data):
    if famcode == "Father" or famcode == "Mother":
        famcode = famcode[0]

    if not path.exists(file_path + "/data2file.csv"):
        with open(file_path + "/data2file.csv", 'w') as f1:
            writer1 = csv.writer(f1)
            writer1.writerow(["FAMCODE", "BIRTHSEQ", "A1", "A2", "B1", "B2", "C1", "C2", "DRB11", "DRB12", "DQB11", "DQB12"])

    with open(file_path + "/data2file.csv", 'a+') as f2:
        writer2 = csv.writer(f2)
        # pay attention the child must be 1, 2,...
        row = ["1", famcode, A1_data, A2_data, B1_data, B2_data, C1_data, C2_data, DRB11_data, DRB12_data, DQB11_data, DQB12_data]
        writer2.writerow(row)


# process basic file the valid format (that the code know how to work with it)
# basic file could arrives in 2 option: 1. data in columns, but not enough (2 columns for each allele, not 4)
# 2. data in gl string (ex: 1,F, A*02+A*03^B*01:01+B*03:01^...)
def basic_csv2format(input_path, output_path):
    wf = open(output_path, 'w')
    writer = csv.writer(wf)
    writer.writerow(
        ['FAMCODE', 'BIRTHSEQ', 'HLA-A*1', 'HLA-A*2', 'A1', 'A2', 'B_STAR1', 'B_STAR2', 'B1', 'B2', 'C_STAR1',
         'C_STAR2', 'C1', 'C2', 'DRB11', 'DRB12', 'DR1', 'DR2', 'DQB11', 'DQB12', 'DQ1', 'DQ2'])
    for line in open(input_path, 'r'):
        csv_row = line.split(',')
        row_to_output = ['']*22
        if csv_row[0].isdigit():
            row_to_output[0] = csv_row[0]
            row_to_output[1] = csv_row[1]

            if "^" in csv_row[2]:  # the input file is in format gl string
                match_als2col = {"A": 2, "B": 6, "C": 10, "DRB1": 14, "DQB1": 18}
                all_als = csv_row[2].split("^")
                for pair_als in all_als:
                    al1, al2 = pair_als.split("+")[0], pair_als.split("+")[1]
                    for key in match_als2col.keys():
                        if key in al1:
                            col_index = match_als2col[key]
                    j = 0
                    for al in [al1, al2]:
                        if "/" in al:  # ambiguity: A*03:01/A*03:02/...
                            row_to_output[col_index + 2 + j] = al.split("*")[1].split(":")[0].replace('\n', '')
                        elif ":" in al:  # high res: A*05:01
                            row_to_output[col_index + j] = al.split("*")[1].replace('\n', '')
                        else:  # low res: A*04
                            row_to_output[col_index + 2 + j] = al.split("*")[1].replace('\n', '')
                        j += 1

            else:  # the input is in format of columns
                for i in range(2, 12):
                    if i % 2 == 0:
                        far_column = i * 2
                    else:
                        far_column = i * 2 - 1
                    if ":" in csv_row[i]:  # high res
                        row_to_output[far_column - 2] = csv_row[i].replace('\n', '')
                    else:  # low res
                        row_to_output[far_column] = csv_row[i].replace('\n', '')
            writer.writerow(row_to_output)


def split_gl(gl_str):
    gl_str = gl_str.strip(' ')
    pairs = gl_str.split("^")
    d = {"A": [], "B": [], "C": [], "DRB1": [], "DQB1": []}
    for pair in pairs:
        # example: pair = A*02:01+A*03:01; allele = A; val1 = 02:01; val2 = 03:01
        allele = pair.split("*")[0]
        val1 = pair.split("+")[0].split("*")[1]
        val2 = pair.split("+")[1].split("*")[1]
        d[allele].extend([val1, val2])
    for key in d.keys():
        if d[key] == []:
            d[key] = ["", ""]
    return d["A"][0], d["B"][0], d["C"][0], d["DRB1"][0], d["DQB1"][0], d["A"][1], d["B"][1], d["C"][1], d["DRB1"][1], d["DQB1"][1]


# create gl string from split data
#  02:01, 30:04, 40:01, 51:22 ... 03:02, 03:02 --->
#  A*02:01+A*30:04^B*40:01+B*51:22^C*03+C*14^DRB1*04:01+DRB1*14:02^DQB1*03:02+DQB1*03:02
def create_glstring(A1, A2, B1, B2, C1, C2, DRB1, DRB2, DQB1, DQB2):
    dict = {A1: "A", A2: "A", B1: "B", B2: "B", C1: "C", C2: "C", DRB1: "DRB1", DRB2: "DRB1", DQB1: "DQB1", DQB2: "DQB1"}
    gl_lst = []
    for pair in [[A1, A2], [B1, B2], [C1, C2], [DRB1, DRB2], [DQB1, DQB2]]:
        pair_lst = []
        for single in pair:
            if single == "":
                continue
            single = str(single)
            full_single = dict[single] + '*' + str(single)  # e.g: A --> A*02:08
            pair_lst.append(full_single)
        pair_gl = '+'.join(pair_lst)    # e.g: [A*02:08, A*03:01] --> A*02:08+A*03:01
        if pair_gl:  # if not empty
            gl_lst.append(pair_gl)
    gl_string = '^'.join(gl_lst)
    return gl_string
















