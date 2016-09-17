# -*- coding: utf-8 -*-
"""
Copyright (C) 2014 TopCoder Inc., All Rights Reserved.

This module is a hack for the MAC OS X dylibs loading issue.
Users should run this script under the "lib" folder of mesh-generator.

@author:  caoweiquan322
@version: 1.0
"""
import shutil
import os

if __name__ == '__main__':
    lines = os.listdir(os.getcwd())
    for line in lines:
        if line.endswith('.5.10.1.dylib'):
            print(line[:-13])
            filename = line[:-13]
            os.remove('./'+filename+'.5.10.dylib')
            os.rename('./'+line, './'+filename+'.5.10.dylib')
