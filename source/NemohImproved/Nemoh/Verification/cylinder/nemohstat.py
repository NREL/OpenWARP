#!/usr/bin/env python
#__author__ = 'yedtoss'

import subprocess
import timeit
import sys,os

preprocessor = ""
solver = ""
num_iterations = 1


def run():
    subprocess.call(preprocessor, shell=True)
    subprocess.call(solver,  shell=True)


if __name__ == "__main__":

    if len(sys.argv) < 3:
        print("Error you need to specify a least 2 arguments")
        print("Usage: ")
        print("python nemohstat.py preProcessorPath SolverPath [num_repeat]")
        sys.exit()

    preprocessor = sys.argv[1]
    solver = sys.argv[2]

    if not os.path.exists(preprocessor) or not os.path.isfile(preprocessor):
        print("Error, preProcessor file does not exists")
        sys.exit()
    if not os.path.exists(solver) or not os.path.isfile(solver):
        print("Error, preProcessor file does not exists")
        sys.exit()
    if len(sys.argv) >= 4:
        num_iterations = int(sys.arg[3])

    wall_time = timeit.timeit(run, number=num_iterations)

    print(" Best Elapsed wall time is {} seconds".format(wall_time))

