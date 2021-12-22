# The purpose of this file is to perform a series of data manipuation and processing commands to particle tracking data in bulk. 
# In particular, functions in this file import particle tracking and ram pressure data, join them as necessary, calculate kinetic 
# and potential energies of particles, classify particles as disk vs. halo, identify ejected or expulsed particles, and more. 
# The reason these functions are written here is so that we can ensure that we are using the same data processing procedures  
# throughout the analysis and not have to repeat this code for each analysis component. 
#
# ____________________________________________________________________________________________
# Code credit to Hollis Akins 2021;
# Github permalink: https://github.com/hollisakins/Justice_League_Code/blob/ 
#                    e049137edcfdc9838ebb3cf0fcaa4ee46e977cec/Analysis/RamPressure/analysis.py
# ____________________________________________________________________________________________
# Last revised: 12 Dec. 2021

import pynbody
import pandas as pd
import numpy as np
import pickle

from base import *



def calc_angles(d):
    '''
    -> Calculates exit angles (angle made between velocity vec of gas particle and velocity vec of host galaxy in satellite
        rest frame) for selected gas particles.
    '''
    #--------------------------------#
    
    # get gas particle velocity
    v = np.array([d.vx,d.vy,d.vz])

    # get velocity of CGM wind (host velocity relative to satellite)
    v_sat = np.array([d.sat_vx,d.sat_vy,d.sat_vz])
    v_host = np.array([d.host_vx,d.host_vy,d.host_vz])
    v_rel = v_host - v_sat # we want the velocity of the host in the satellite rest frame

    # take the dot product and get the angle, in degrees
    v_hat = v / np.linalg.norm(v)
    v_rel_hat = v_rel / np.linalg.norm(v_rel)
    angle = np.arccos(np.dot(v_hat,v_rel_hat)) * 180/np.pi

    d['angle'] = angle
        
    return d



def calc_angles_tidal(d):
    # get gas particle velocity
    v = np.array([d.vx,d.vy,d.vz])

    # instead of CGM velocity, get vector pointing from satellite to host (i.e. host position in the satellite rest frame) 
    r_sat =np.array([d.sat_Xc,d.sat_Yc,d.sat_Zc])
    r_host = np.array([d.host_Xc,d.host_Yc,d.host_Zc])
    r_rel = r_host - r_sat

    # take the dot product and get the angle, in degrees
    v_hat = v / np.linalg.norm(v)
    r_rel_hat = r_rel / np.linalg.norm(r_rel)
    angle = np.arccos(np.dot(v_hat,r_rel_hat)) * 180/np.pi

    d['angle_tidal'] = angle
        
    return 



def vec_to_xform(vec):
    '''
    -> Aligns Pynbody coordinate system to a specified vector 'vec'.
    '''
    #--------------------------------#
    
    vec_in = np.asarray(vec)
    vec_in = vec_in / np.sum(vec_in ** 2).sum() ** 0.5
    vec_p1 = np.cross([1, 0, 0], vec_in)
    vec_p1 = vec_p1 / np.sum(vec_p1 ** 2).sum() ** 0.5
    vec_p2 = np.cross(vec_in, vec_p1)
    matr = np.concatenate((vec_p2, vec_in, vec_p1)).reshape((3, 3))
    return matr



