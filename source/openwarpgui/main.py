# -*- coding: utf-8 -*-
"""
Copyright (C) 2014-2015 TopCoder Inc., All Rights Reserved.

This module defined controllers for CherryPy.

Updated since version 1.1:
    1. Updated the method of launching browser to  be based on lamda. This is more brief.
    2. Updated the CherryPy gloabal configuration so that session is supported.

Updated since version 1.2: Merge Code and Update GUI
    1. Integrate New Nemoh using hdf5 and python.

@author:  caoweiquan322, TCSASSEMBLER
@version: 1.2
"""
from openwarp.settings import *

import webbrowser
import cherrypy
import time
import threading
import subprocess

# Entry-point of the whole application
if __name__ == '__main__':
    # Compile python module if it was not compiled
    subprocess.call(['python', 'setup.py', 'build_ext', '--inplace'], cwd='nemoh')
    # Start browser starter thread
    # The browser will start after 2 seconds
    threading.Timer(2,
                    lambda: webbrowser.open('http://127.0.0.1:' +
                                            str(WEB_SERVER_PORT) +
                                            '/index.html')).start()
    
    # Start up web server
    from openwarp.web import WebController
    cherrypy.quickstart(WebController(), '', {
        'global': {
            'server.socket_host': '127.0.0.1',
            'server.socket_port': WEB_SERVER_PORT,
            'tools.sessions.on': True,
            'tools.sessions.timeout': 60000000,
            #'server.socket_timeout': 60000000
        },
        '/': {
            'tools.staticdir.dir' : os.path.join(os.path.dirname(os.path.abspath(__file__)), 'static'),
            'tools.staticdir.on' : True,
            'response.timeout': 60000000
        }
    })
