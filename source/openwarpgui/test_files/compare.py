#!/usr/bin/env python
#__author__ = 'yedtoss'

import os.path
from operator import sub,abs
import sys
from os import walk
from math import isnan
import numpy as np


def isfloat(str1):
    """
    Check if a string is a float
    """
    try:
        float(str1)
    except ValueError:
        return False
    return True


def compute(filename1, filename2):
    """
    Compute the MAE between real numbers in two files.
    It drops all non number. It use 0 for NaN
    """

    x1 = []
    x2 = []

    if os.path.exists(filename1) and os.path.isfile(filename1):
        f = open(filename1, 'r')
        for line in f:
                for word in line.split():
                    if isfloat(word):
                        x1.append(0 if isnan(float(word)) else float(word))
        f.close()

    if os.path.exists(filename2) and os.path.isfile(filename2):
        f = open(filename2, 'r')
        for line in f:
                for word in line.split():
                    if isfloat(word):
                        x2.append(0 if isnan(float(word)) else float(word))
        f.close()

    m = min(len(x1), len(x2))
    score = sum(map(abs, map(sub, x1[:m], x2[:m])))
    print(len(x1))
    print(len(x2))
    if m > 0:
        score = score/m

    print("MAE between {} and {} is {}".format(filename1, filename2, score))

    #np.savetxt(filename1+ '_correct', x1)
    #np.savetxt(filename2+ '_correct', x2)

    return score


if __name__ == "__main__":

    if len(sys.argv) < 3:
        print("Error you need to specify a least two arguments")
        print("Usage: ")
        print("python compare.py originalfile newfile")
        print("python compare.py originaldirectory newdirectory [basestring]")
        sys.exit()

    name1 = sys.argv[1]
    name2 = sys.argv[2]
    base = ""
    if len(sys.argv) >= 4:
        base = sys.argv[3]
    if not os.path.exists(name1):
        print("Error: {} does not exist".format(name1))
        sys.exit()
    if not os.path.exists(name2):
        print("Error: {} does not exist".format(name2))
        sys.exit()

    if os.path.isfile(name1) and os.path.isdir(name2):
        print("We can't compare file and directory")
        sys.exit()

    if os.path.isfile(name2) and os.path.isdir(name1):
        print("We can't compare file and directory")
        sys.exit()

    if os.path.isfile(name1) and os.path.isfile(name2):
        compute(name1, name2)

    if os.path.isdir(name1) and os.path.isdir(name2):

        list1 = []
        list2 = []
        root1 = ""
        root2 = ""
        total = 0.0
        num = 0

        for (root, dirs, files) in walk(name1):
            list1.extend(files)
            root1 = root
            break

        for (root, dirs, files) in walk(name2):
            list2.extend(files)
            root2 = root
            break

        for name in list1:
            if name in list2:

                if name.startswith(base):
                    total += compute(os.path.join(root1, name), os.path.join(root2, name))
                    num += 1

        if num > 0:
            total = total/num

        print("Total Average MAE between {} and {} is {}".format(name1, name2, total))
