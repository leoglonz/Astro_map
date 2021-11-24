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



def calc_adv_accreted(sim, haloid, save=True, verbose=True):
    ''' 
    Computing accreted gas particles using method 1: a particle is recorded as accreted 
    if it was present in the 'discharged' dataset prior to the accretion event.
    '''
    
    import tqdm

    discharged, accreted = read_discharged(sim, haloid, verbose=verbose)  

    if verbose: print(f'Now computing adv. accreted particles for {sim}-{haloid}...')
    adv_accreted = pd.DataFrame() # gas accreted following a discharge event.

    pids = np.unique(accreted.pid)



    for pid in tqdm.tqdm(pids):
        dis = discharged[discharged.pid==pid]
        acc = accreted[accreted.pid==pid]

        dTime = np.asarray(dis.time)
        aTime = np.asarray(acc.time)


        # Case 1: removing initial accretion event if it does not correspond to a discharge; else no alterations.
        # check to ensure that our discharge event actually has an accretion.
        if (len(aTime) != 0) and (len(dTime) != 0):
            if (aTime[0] < dTime[0]):
                aCache = acc[1:]

            else:
                aCache = acc
            

        adv_accreted = pd.concat([adv_accreted, aCache])
        
        
    if save:
        key = f'{sim}_{str(int(haloid))}'
        filepath = '/home/lonzaric/astro_research/Stellar_Feedback_Code/SNeData/adv_accreted.hdf5'
        print(f'Saving {key} adv. accreted particle dataset to {filepath}')
        adv_accreted.to_hdf(filepath, key=key)
        
        
        
    print(f'> Returning (adv. accreted) datasets <')

    return adv_accreted
            
    
    
def read_accreted():
    adv_accreted = pd.DataFrame()

    keys = get_keys()

    for i,key in enumerate(keys):
        i += 1
        sim = key[:4]
        haloid = int(key[5:])
        adv_accreted1 = pd.read_hdf('/home/lonzaric/astro_research/Stellar_Feedback_Code/SNeData/adv_accreted.hdf5',\
             key=key)
        adv_accreted1['key'] = key
        adv_accreted = pd.concat([adv_accreted, adv_accreted1])


    print(f'> Returning (adv. accreted) for all available satellites <')
    return adv_accreted

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
 