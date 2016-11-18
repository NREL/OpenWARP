# -*- coding: utf-8 -*-
"""
This is module that represents the helper methods.

Updated since version 1.1:
    1. Added check_path_exist(), check_is_directory() and check_is_file().

Updated since version 1.2 (OpenWarp - Add Logging Functionality) :
    Added support for logging
"""

__author__ = "caoweiquan322, yedtoss"
__copyright__ = "Copyright (C) 2014-2016 TopCoder Inc. All rights reserved."
__version__ = "1.2"

import traceback
import os

def check_not_none_nor_empty(val, name):
    '''
    Check if the given value is None or empty string.

    @param val: the given value to check
    @param name: name of val
    @raise TypeError: if val is not of type string
    @raise ValueError: if val is None or empty string
    '''
    if val == None:
        raise ValueError('Object ' + name + ' should not be None.')
    if not isinstance(val, str) and not isinstance(val, unicode):
        raise TypeError('Object ' + name + ' should be a string.')
    if len(val.strip()) == 0:
        raise ValueError('Object ' + name + ' should not empty string.')

def check_path_exist(val, name):
    '''
    Check if the given value is a legal file or directory path.

    @param val: the given value to check
    @param name: name of val
    @raise ValueError: if val is not a legal file or directory path
    '''
    if not os.path.exists(val):
        raise ValueError('Path "' + val + '" does not exist.')

def check_is_directory(val, name):
    '''
    Check if the given value is a legal directory path.

    @param val: the given value to check
    @param name: name of val
    @raise ValueError: if val is not a legal directory path
    '''
    check_path_exist(val, name)
    if not os.path.isdir(val):
        raise ValueError('Path "' + val + '" is not a legal directory.')

def check_is_file(val, name):
    '''
    Check if the given value is a legal file path.

    @param val: the given value to check
    @param name: name of val
    @raise ValueError: if val is not a legal file path
    '''
    check_path_exist(val, name)
    if not os.path.isfile(val):
        raise ValueError('Path "' + val + '" is not a legal file.')

def check_type_value(val, name, expected_type, allow_none):
    '''
    Check if the given value is of expected type. And also check if the val is None.

    @param val: the given value to check
    @param name: name of val
    @param expected_type: the expected type
    @param allow_none: whether the val is allowed to be None
    @raise TypeError: if val is not of expected type
    @raise ValueError: if val is None while not allow None
    '''
    if val == None and not allow_none:
        raise ValueError('Object ' + name + ' should not be None.')
    if not isinstance(val, expected_type):
        raise TypeError('Object ' + name + ' should be of ' + str(expected_type) + '.')

def log_entrance(logger, signature, parasMap):
    '''
    Logs for entrance into public methods at DEBUG level.

    @param logger: the logger object
    @param signature: the method signature
    @param parasMap: the passed parameters
    '''
    logger.debug('[Entering method ' + signature + ']')
    if parasMap != None and len(parasMap.items()) > 0:
        paraStr = '[Input parameters['
        for (k,v) in parasMap.items():
            paraStr += (str(k) + ':' + str(v) + ', ')
        paraStr += ']]'
        logger.debug(paraStr)

def log_exit(logger, signature, parasList):
    '''
    Logs for exit from public methods at DEBUG level.

    @param logger: the logger object
    @param signature: the method signature
    @param parasList: the objects to return
    '''
    logger.debug('[Exiting method ' + signature + ']')
    if parasList != None and len(parasList) > 0:
        logger.debug('[Output parameter ' + str(parasList) + ']')

def log_exception(logger, signature, e):
    '''
    Logging exception at ERROR level.

    @param logger: the logger object
    @param signature: the method signature
    @param e: the error
    '''
    # This will log the traceback.
    logger.error('[Error in method ' + signature + ': Details ' + str(e) + ']')
    logger.error(' Error stack:')
    logger.error(traceback.format_exc())
    return e
