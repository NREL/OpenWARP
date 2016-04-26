# -*- coding: utf-8 -*-
"""
This module defined controllers for CherryPy.

Updated since version 1.1:
    1. Added support for postprocess and visualization.
    2. Added support for quit.

Updated since version 1.2 (OpenWarp - Add Logging Functionality):
    Added support for logging
"""

__author__ = "caoweiquan322, yedtoss"
__copyright__ = "Copyright (C) 2014-2016 TopCoder Inc. All rights reserved."
__version__ = "1.2"

import cherrypy
from openwarp import services
from openwarp.services import *
import os
import threading
import json
import logging
import helper

class WebController:
    '''
    This class exposes HTTP services for the frontend HTML to consume using AJAX.
    '''

    def __init__(self):
        self.logger = logging.getLogger(__name__ + '.WebController')

    @cherrypy.expose
    @cherrypy.tools.json_out()
    def apply_configuration(self, **kwargs):
        '''
        Apply the application wide configuration.

        @param self: the class instance itself
        @param kwargs: the other arguments
        @return: the response as a dictionary, will be serialized to JSON by CherryPy.
        '''
        signature = __name__ + '.WebController.apply_configuration()'
        helper.log_entrance(self.logger, signature, kwargs)

        try:

            # Prepare meshing directory
            self.logger.info('Applying configuration ...')

            # Call generate_mesh service
            ret = {
                'log': services.apply_configuration(ConfigurationParameters(**kwargs))
            }
            helper.log_exit(self.logger, signature, [ret])
            return ret

        except (TypeError, ValueError) as e:
            helper.log_exception(self.logger, signature, e)
            # Error with input, respond with 400
            cherrypy.response.status = 400
            ret = { 'error' : str(e) }
            helper.log_exit(self.logger, signature, [ret])
            return ret
        except Exception as e:
            helper.log_exception(self.logger, signature, e)
            # Server internal error, respond with 500
            cherrypy.response.status = 500
            ret = { 'error' : str(e) }
            helper.log_exit(self.logger, signature, [ret])
            return ret


    @cherrypy.expose
    @cherrypy.tools.json_out()
    def generate_mesh(self, **kwargs):
        '''
        Launch Mesh Generator to generate mesh.

        @param self: the class instance itself
        @param kwargs: the other arguments
        @return: the response as a dictionary, will be serialized to JSON by CherryPy.
        '''
        signature = __name__ + '.WebController.generate_mesh()'
        helper.log_entrance(self.logger, signature, kwargs)
        
        try:
            # Prepare meshing directory
            self.logger.info('Preparing meshing directory')
            meshing_dir = services.prepare_dir('meshing_')
            self.logger.info('Meshing files will be located at ' + str(meshing_dir))

            # Call generate_mesh service
            ret = {
                'log' : services.generate_mesh(meshing_dir, MeshingParameters(**kwargs))
            }
            helper.log_exit(self.logger, signature, [ret])
            return ret
        except (TypeError, ValueError) as e:
            helper.log_exception(self.logger, signature, e)
            # Error with input, respond with 400
            cherrypy.response.status = 400
            ret = { 'error' : str(e) }
            helper.log_exit(self.logger, signature, [ret])
            return ret
        except Exception as e:
            helper.log_exception(self.logger, signature, e)
            # Server internal error, respond with 500
            cherrypy.response.status = 500
            ret = { 'error' : str(e) }
            helper.log_exit(self.logger, signature, [ret])
            return ret

    @cherrypy.expose
    @cherrypy.tools.json_out()
    def simulate(self, json_str):
        '''
        Run simulation.

        @param self: the class instance itself
        @param json_str: the json string posted by client
        @return: the response as a dictionary, will be serialized to JSON by CherryPy.
        '''
        signature = __name__ + '.WebController.simulate()'
        helper.log_entrance(self.logger, signature, {'json_str' : json_str})
        
        try:
            # Prepare simulation directory
            simulation_dir = services.prepare_dir('simulation_')
            cherrypy.session['simulation_dir'] = simulation_dir
            cherrypy.session['simulation_done'] = False
            # Call simulate service
            ret = {
                'log' : services.simulate(simulation_dir, self.construct_simulation_parameters(json_str))
            }
            cherrypy.session['simulation_done'] = True
            # Set postprocess flag to False if a new simulation has been done successfully.
            cherrypy.session['postprocess_done'] = False
            helper.log_exit(self.logger, signature, [ret])
            return ret
        except (TypeError, ValueError) as e:
            helper.log_exception(self.logger, signature, e)
            # Error with input, respond with 400
            cherrypy.response.status = 400
            ret = { 'error' : str(e) }
            helper.log_exit(self.logger, signature, [ret])
            return ret
        except Exception as e:
            helper.log_exception(self.logger, signature, e)
            # Server internal error, respond with 500
            cherrypy.response.status = 500
            ret = { 'error' : str(e) }
            helper.log_exit(self.logger, signature, [ret])
            return ret

    @cherrypy.expose
    @cherrypy.tools.json_out()
    def upload_file(self, uploadedFile):
        '''
        Upload a file via AJAX request, the file will be created in temporary directory and the full path will
        be sent back as JSON response.

        @param self: the class instance itself
        @param uploadedFile: the uploaded file
        @return: the response as a dictionary, will be serialized to JSON by CherryPy.
        '''
        signature = __name__ + '.WebController.upload_file()'
        helper.log_entrance(self.logger, signature, {'uploadedFile' : uploadedFile})
        
        try:
            temp_dir = os.path.join(TEMP_DATA_DIRECTORY, uuid.uuid1().hex)
            os.mkdir(temp_dir)
            filepath = os.path.join(temp_dir, uploadedFile.filename)
            
            # We must use 'wb' mode here in case the uploaded file is not ascii format.
            with open(filepath, 'wb') as output:
                while True:
                    data = uploadedFile.file.read(1024)
                    if data:
                        output.write(data)
                    else:
                        break
            try:
                with open(filepath, 'r') as input:
                    points, panels = self.determine_points_panels(input)
                    ret = {
                        'filepath' : filepath,
                        'points' : points,
                        'panels' : panels
                    }
                    helper.log_exit(self.logger, signature, [ret])
                    return ret
            except Exception as e:
                helper.log_exception(self.logger, signature, e)
                ret = { 'filepath' : filepath }
                helper.log_exit(self.logger, signature, [ret])
                return ret
        except Exception as e:
            helper.log_exception(self.logger, signature, e)
            # Server internal error, respond with 500
            cherrypy.response.status = 500
            ret = { 'error' : str(e) }
            helper.log_exit(self.logger, signature, [ret])
            return ret

    @cherrypy.expose
    @cherrypy.tools.json_out()
    def postprocess(self, json_str):
        '''
        Run post-processing.

        @param self: the class instance itself
        @param json_str: the json string posted by client
        @return: the response as a dictionary, will be serialized to JSON by CherryPy.
        '''
        signature = __name__ + '.WebController.postprocess()'
        helper.log_entrance(self.logger, signature, {'json_str': json_str})
        # Set session variable postprocess_done to False by default.
        cherrypy.session['postprocess_done'] = False
        
        try:
            if not cherrypy.session.has_key('simulation_done') or not cherrypy.session['simulation_done']:
                # simulation must be run first
                cherrypy.response.status = 400
                ret = { 'error' : 'Simulation must be run first.' }
                helper.log_exit(self.logger, signature, [ret])
                return ret
            else:
                # Call post-processing service
                ret = {
                    'log' : services.postprocess(cherrypy.session['simulation_dir'],
                                                 self.construct_postprocess_parameters(json_str))
                }
                cherrypy.session['postprocess_done'] = True
                helper.log_exit(self.logger, signature, [ret])
                return ret
        except (TypeError, ValueError) as e:
            helper.log_exception(self.logger, signature, e)
            # Error with input, respond with 400
            cherrypy.response.status = 400
            ret = { 'error' : str(e) }
            helper.log_exit(self.logger, signature, [ret])
            return ret
        except Exception as e:
            helper.log_exception(self.logger, signature, e)
            # Server internal error, respond with 500
            cherrypy.response.status = 500
            ret = { 'error' : str(e) }
            helper.log_exit(self.logger, signature, [ret])
            return ret

    @cherrypy.expose
    @cherrypy.tools.json_out()
    def visualize(self):
        '''
        Launch ParaView to visualize simulation results.

        @param self: the class instance itself
        @return: the response as a dictionary, will be serialized to JSON by CherryPy.
        '''
        signature = __name__ + '.WebController.visualize()'
        helper.log_entrance(self.logger, signature, None)
        
        try:
            if not cherrypy.session.has_key('simulation_done') or not cherrypy.session['simulation_done']:
                # simulation must be run first
                cherrypy.response.status = 400
                ret = { 'error' : 'Simulation must be run first.' }
                helper.log_exit(self.logger, signature, [ret])
                return ret
            elif not cherrypy.session.has_key('postprocess_done') or not cherrypy.session['postprocess_done']:
                # postprocess must be run first
                cherrypy.response.status = 400
                ret = { 'error' : '"SAVE AS TECPLOT" must be run right after a successful simulation.' }
                helper.log_exit(self.logger, signature, [ret])
                return ret
            else:
                # Call visualize service
                services.visualize(cherrypy.session['simulation_dir'])
                helper.log_exit(self.logger, signature, None)
                return {}
        except (TypeError, ValueError) as e:
            helper.log_exception(self.logger, signature, e)
            # Error with input, respond with 400
            cherrypy.response.status = 400
            ret = { 'error' : str(e) }
            helper.log_exit(self.logger, signature, [ret])
            return ret
        except Exception as e:
            helper.log_exception(self.logger, signature, e)
            # Server internal error, respond with 500
            cherrypy.response.status = 500
            ret = { 'error' : str(e) }
            helper.log_exit(self.logger, signature, [ret])
            return ret

    @cherrypy.expose
    @cherrypy.tools.json_out()
    def quit(self):
        '''
        Quit the application by shutting down the CherryPy server.

        @param self: the class instance itself
        @return: the response as a dictionary, will be serialized to JSON by CherryPy.
        '''
        signature = __name__ + '.WebController.quit()'
        helper.log_entrance(self.logger, signature, None)
        
        # Quit after sending response
        threading.Timer(2, lambda: os._exit(0)).start()
        helper.log_exit(self.logger, signature, None)
        return {}
    
    def construct_simulation_parameters(self, json_str):
        '''
        Construct the simulation parameters from json string.

        @param self: the class instance itself
        @param json_str: the json string to parse
        @return: the parsed SimulationParameters object
        '''
        # Since this is a internal method. The parameters won't be logged.
        json_obj = json.JSONDecoder().decode(json_str)
        para = SimulationParameters(**json_obj)
        if para.floating_bodies is not None:
            new_bodies = []
            for body in para.floating_bodies:
                new_bodies.append(FloatingBody(**body))
            del para.floating_bodies[:]
            para.floating_bodies.extend(new_bodies)
        return para

    def construct_postprocess_parameters(self, json_str):
        # Since this is a internal method. The parameters won't be logged.
        json_obj = json.JSONDecoder().decode(json_str)
        para = PostprocessingParameters(**json_obj)
        return para

    def determine_points_panels(self, dat_file):
        '''
        Determines the number of points and panels of a mesh file.

        @param self: the class instance itself
        @param dat_file: the mesh file to parse
        @return: the number of points and panels of a mesh file
        @raise Exception: if the file is not expected format
        '''
        # Since this is a internal method. The parameters won't be logged.
        lines = dat_file.readlines()
        num_lines = 0
        zero_line1 = 0
        zero_line2 = 0
        succeed = False
        for line in lines[1:]:
            if len(line.strip()) > 0:
                zero_line1 = zero_line1 + 1
                if line.strip().startswith('0'):
                    succeed = True
                    break
        if not succeed:
            raise Exception('Zero line 1 not found.')
        for line in lines[(zero_line1 + 1):]:
            if len(line.strip()) > 0:
                zero_line2 = zero_line2 + 1
                if line.strip().startswith('0'):
                    succeed = True
                    break
        if not succeed:
            raise Exception('Zero line 2 not found.')
        return zero_line1 - 1, zero_line2 - 1
            
