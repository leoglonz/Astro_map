# The purpose of this file is to perform a series of data manipuation and processing commands to particle tracking data in bulk. 
# In particular, functions in this file import particle tracking and sne-heating data; calculate necessary particle attributes;
# classify particles as disk vs. halo; identify predischarged, discharged, and acrreted particles; etc.
# We use this . 
import pynbody
import pandas as pd
import numpy as np
import pickle
from base import *

# from analysis import *


#####################################################################
### Used to figure out sufficient way to compute accreted dataset ###
#####################################################################



def get_keys():
    path1 = '/home/lonzaric/astro_research/Stellar_Feedback_Code/SNeData/ejected_particles.hdf5'
    with pd.HDFStore(path1) as hdf:
        keys = [k[1:] for k in hdf.keys()]
    print(*keys)
    return keys



def read_tracked_particles(sim, haloid, verbose=False):
    
    if verbose: print(f'Loading tracked particles for {sim}-{haloid}...')
    
    key = f'{sim}_{str(int(haloid))}'

    # import the tracked particles dataset
    path1 = '/home/lonzaric/astro_research/Stellar_Feedback_Code/SNeData/tracked_particles_v2.hdf5'
    data = pd.read_hdf(path1, key=key)
    
    time = np.unique(data.time)
    dt = time[1:]-time[:-1]
    dt = np.append(dt[0], dt)
    dt = dt[np.unique(data.time, return_inverse=True)[1]]
    data['dt'] = dt
    
    
    if verbose: print('Successfully loaded')
    
    r_gal = np.array([])
    for t in np.unique(data.time):
        d = data[data.time==t]
        r_gas = np.mean(d.sat_r_gas)
        r_half = np.mean(d.sat_r_half)
        rg = np.max([r_gas,r_half])

        if np.isnan(rg):
            rg = r_gal_prev

        if verbose: print(f't = {t:1f} Gyr, satellite R_gal = {rg:.2f} kpc')
        r_gal = np.append(r_gal,[rg]*len(d))

        r_gal_prev = rg

    data['r_gal'] = r_gal
    
    r_gal_prev = 0
    r_gal = np.array([])
    for t in np.unique(data.time):
        d = data[data.time==t]
        r_gas = np.mean(d.host_r_gas)
        r_half = np.mean(d.host_r_half)
        rg = np.max([r_gas,r_half])

        if np.isnan(rg):
            rg = r_gal_prev

        if verbose: print(f't = {t:1f} Gyr, host R_gal = {rg:.2f} kpc')
        r_gal = np.append(r_gal,[rg]*len(d))

        r_gal_prev = rg

    data['host_r_gal'] = r_gal
    
    thermo_disk = (np.array(data.temp) < 1.2e4) & (np.array(data.rho) > 0.1)
    
    in_sat = np.array(data.in_sat)
    other_sat = np.array(data.in_other_sat)
    in_host = np.array(data.in_host) & ~in_sat & ~other_sat
    
    sat_disk = in_sat & thermo_disk
    sat_halo = in_sat & ~thermo_disk
    
    host_disk = in_host & thermo_disk
    host_halo = in_host & ~thermo_disk
    
    IGM = np.array(data.in_IGM)
    
    # basic location classifications.
    data['sat_disk'] = sat_disk
    data['sat_halo'] = sat_halo
    data['host_disk'] = host_disk
    data['host_halo'] = host_halo
    data['other_sat'] = other_sat
    data['IGM'] = IGM

    return data



def read_discharged(sim, haloid, verbose=False):
    discharged = pd.DataFrame()
    
    key = f'{sim}_{str(int(haloid))}'
    
    path = '/home/lonzaric/astro_research/Stellar_Feedback_Code/SNeData/discharged_particles.hdf5'
    discharged = pd.read_hdf(path, key=key)
    
    path = '/home/lonzaric/astro_research/Stellar_Feedback_Code/SNeData/dsrg_accreted_particles.hdf5'
    accreted = pd.read_hdf(path, key=key)
    
    return discharged, accreted



def calc_accreted(sim, haloid, save=True, verbose=True):
    ''' 
    Computing accreted gas particles using method 1: a particle is recorded as accreted 
    if it was present in the 'discharged' dataset prior to the accretion event.
    '''
    
    import tqdm
    
    tracked = read_tracked_particles(sim, haloid, verbose=verbose)# all gas at all recorded times.
    discharged = read_discharged(sim, haloid, verbose=verbose) # all gas ejected from the disk (can include repeat events).
    
#     if verbose: print(f'Now computing discharged particles for {sim}-{haloid}...')
    dsrg_accreted = pd.DataFrame() # gas accreted following a discharge event.
    
    
    # picking out all unique particles that were tracked in a the dwarf galaxy.
    pids = np.unique(tracked.pid) 
    dsrg_pids = np.unique(discharged.pid)
    
    
    
    
    for pid in tqdm.tqdm(dsrg_pids): # iterating over each unique pid.
        data = tracked[tracked.pid==pid] # picking out all instances of the particular particle (each has same # of timesteps).
        dsrg_data = discharged[discharged.pid==pid]
    
        
        sat_disk = np.array(data.sat_disk, dtype=bool)
        in_sat = np.array(data.in_sat, dtype=bool)
        outside_disk = ~sat_disk
        
        time = np.array(data.time, dtype=float)
        
                
            
         
            
            
        for i,t2 in enumerate(time[1:]): # iterating over all recorded timesteps for the particle.
                i += 1
                    
                if ( and sat_disk[i]):
                    _in_ = data[time==t2].copy()
                    dsrg_accreted = pd.concat([dsrg_accreted, _in_])
                 

                
                
                
                
                
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