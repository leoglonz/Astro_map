# The purpose of this file is to perform a series of data manipuation and processing commands to particle tracking data in bulk. 
# In particular, functions in this file import particle tracking and sne-heating data; calculate necessary particle attributes;
# classify particles as disk vs. halo; identify predischarged, discharged, and acrreted particles; etc.
# We use this . 
import pynbody
import pandas as pd
import numpy as np
import pickle
from base import *

from analysis import *


#####################################################################
### Used to figure out sufficient way to compute accreted dataset ###
#####################################################################


def calc_accreted(sim, haloid, save=True, verbose=True):
    import tqdm
    data = read_tracked_particles(sim, haloid, verbose=verbose)

    if verbose: print(f'Now computing discharged particles for {sim}-{haloid}...')
    dsrg_accreted = pd.DataFrame() # gas accreted following a discharge event.

    
    pids = np.unique(data.pid)
    for pid in tqdm.tqdm(pids):
        dat = data[data.pid==pid]

        sat_disk = np.array(dat.sat_disk, dtype=bool)
        in_sat = np.array(data.in_sat, dtype=bool)
        outside_disk = ~sat_disk
        
        time = np.array(dat.time, dtype=float)

        for i,t2 in enumerate(time[1:]):
                i += 1
                
                # stopping condition to avoid enumeration overflow:
                if i == len(time)-1:
                    break 
                    
#                 if sat_disk[i-1] and outside_disk[i] and sat_disk[i+1]:
#                     in_ = dat[time==time[i-1]].copy()
#                     out = dat[time==t2].copy()
#                     predischarged = pd.concat([predischarged, in_])
#                     discharged = pd.concat([discharged, out])
                    
                if (sat_disk[i-1] and outside_disk[i] and sat_disk[i+1]):
                    out = dat[time==t2].copy()
                    dsrg_accreted = pd.concat([dsrg_accreted, out])
                 

    # apply the calc_angles function along the rows of discharged particles.
   
    
    if save:
        filepath = '/home/lonzaric/astro_research/Stellar_Feedback_Code/SNeData/accreted_test.hdf5'
        print(f'Saving {key} accreted particle dataset to {filepath}')
        dsrg_accreted.to_hdf(filepath, key=key)
        
        
    print(f'> Returning (dsrg_accreted) datasets <')

    return dsrg_accreted






keys = ['h148_13','h148_28','h148_37','h148_45','h148_68','h148_80','h148_283',
        'h148_278','h148_329','h229_20','h229_22','h229_23','h229_27','h229_55',
        'h242_24','h242_41','h242_80','h329_33','h329_137']

print('Compiling sim gas into sets (accreted) for the following keys:', keys)

for key in keys:
    sim = str(key[:4])
    haloid = int(key[5:])
    # note that heated is automatically concatenated without double counting, irrespective 
    # of how many times this code is run.
    accreted = calc_accreted(sim, haloid, save=True, verbose=False)