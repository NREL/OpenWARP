# -*- coding: utf-8 -*-
"""
Copyright (C) 2014-2016 TopCoder Inc., All Rights Reserved.

This module defined controllers for CherryPy.

Updated since version 1.1:
    1. Updated the method of launching browser to  be based on lamda. This is more brief.
    2. Updated the CherryPy gloabal configuration so that session is supported.

Updated since version 1.2: Merge Code and Update GUI
    1. Integrate New Nemoh using hdf5 and python.

Changes in version 1.3 (OpenWarp - Add Logging Functionality)
       Added support for logging.

@author:  caoweiquan322, yedtoss
@version: 1.3
"""
from openwarp.settings import *

import webbrowser
import cherrypy
import time
import threading
import subprocess
from nemoh import utility
from nemoh import settings

# Entry-point of the whole application
if __name__ == '__main__':
    utility.setup_logging(default_conf_path=settings.LOGGING_CONFIGURATION_FILE, logging_path=LOG_FILE)
    is_log_dir_exists = os.path.exists(LOG_DIR)
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
    # Stop cherrypy from propagating logs to our log file
    cherrypy.log.error_log.propagate = False
    cherrypy.log.access_log.propagate = False
    cherrypy.quickstart(WebController(), '', {
        'global': {
            'server.socket_host': '127.0.0.1',
            'server.socket_port': WEB_SERVER_PORT,
            'tools.sessions.on': True,
            'tools.sessions.timeout': 60000000,
            'log.access_file': os.path.join(LOG_DIR, 'access.log') if is_log_dir_exists else './access.log',
            'log.error_file': os.path.join(LOG_DIR, 'error.log') if is_log_dir_exists else './error.log',
            'log.screen': False
            #'server.socket_timeout': 60000000
        },
        '/': {
            'tools.staticdir.dir' : os.path.join(os.path.dirname(os.path.abspath(__file__)), 'static'),
            'tools.staticdir.on' : True,
            'response.timeout': 60000000
        }
    })
