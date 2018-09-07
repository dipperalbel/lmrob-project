#!/usr/bin/python

import csv
import sys
import math
with open(sys.argv[1], newline='') as f:
    with open(sys.argv[2], newline='') as f2:
        reader = csv.reader(f, delimiter=';', quoting=csv.QUOTE_NONE)

        reader2 = csv.reader(f2, delimiter=';', quoting=csv.QUOTE_NONE)
        for row, row2 in zip(reader, reader2):
            if '"x"' in row[0]:
                continue
            line = "| {:>14.8f} | {:>19.8f} | {:>10.5e} |".format(float(row[0]), float(row2[0]), math.fabs(float(row2[0]) - float(row[0])))
            print(line)