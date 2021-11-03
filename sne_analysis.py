# The purpose of this file is to perform a series of data manipuation and processing commands to particle tracking data in bulk. 
# In particular, functions in this file import particle tracking and ram pressure data, join them as necessary, calculate kinetic 
# and potential energies of particles, classify particles as disk vs. halo, identify ejected or expulsed particles, and more. 
# The reason these functions are written here is so that we can ensure that we are using the same data processing procedures throughout 
# the analysis and not have to repeat this code for each analysis component. 
import pynbody
import pandas as pd
import numpy as np
import pickle
from base import *