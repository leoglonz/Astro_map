def norm(x):
    '''
    Overview:
        Takes the rudementary norm of the rows of arrays. Quicker than computation via linalg.norm.
    
    Inputs:
    -> 'x', array-like.
    
    Output:
    -> 1xn array with each column entry representing the norm of the row of equivalent index.
    '''
    #-----------------------------#
    
    return np.sum(np.abs(x)**2,axis=-1)**(1./2)



# ------------------------------------------- #

def calc_ejected_expelled(sim, haloid, save=True, verbose=True):
    import tqdm
    data = read_tracked_particles(sim, haloid, verbose=verbose)

    if verbose: print(f'Now computing ejected/expelled particles for {sim}-{haloid}...')
    ejected = pd.DataFrame()
    cooled = pd.DataFrame()
    expelled = pd.DataFrame()
    accreted = pd.DataFrame()
    
    pids = np.unique(data.pid)
    for pid in tqdm.tqdm(pids):
        dat = data[data.pid==pid]

        sat_disk = np.array(dat.sat_disk, dtype=bool)
        sat_halo = np.array(dat.sat_halo, dtype=bool)
        in_sat = np.array(dat.in_sat, dtype=bool)
        outside_sat = ~in_sat

        host_halo = np.array(dat.host_halo, dtype=bool)
        host_disk = np.array(dat.host_disk, dtype=bool)
        IGM = np.array(dat.IGM, dtype=bool)
        other_sat = np.array(dat.other_sat, dtype=bool)
        
        time = np.array(dat.time,dtype=float)



        for i,t2 in enumerate(time[1:]):
                i += 1
                if sat_disk[i-1] and sat_halo[i]:
                    out = dat[time==t2].copy()
                    ejected = pd.concat([ejected, out])
                    
                if sat_halo[i-1] and sat_disk[i]:
                    out = dat[time==t2].copy()
                    cooled = pd.concat([cooled, out])
                    
                if in_sat[i-1] and outside_sat[i]:
                    out = dat[time==t2].copy()
                    if sat_halo[i-1]:
                        out['state1'] = 'sat_halo'
                    elif sat_disk[i-1]:
                        out['state1'] = 'sat_disk'
                        
                    expelled = pd.concat([expelled, out])
                    
                if outside_sat[i-1] and in_sat[i]:
                    out = dat[time==t2].copy()
                    if sat_halo[i]:
                        out['state2'] = 'sat_halo'
                    elif sat_disk[i]:
                        out['state2'] = 'sat_disk'
                        
                    accreted = pd.concat([accreted, out])

    # apply the calc_angles function along the rows of ejected and expelled
    print('Calculating ejection angles')
    ejected = ejected.apply(calc_angles, axis=1)
    print('Calculating expulsion angles')
    expelled = expelled.apply(calc_angles, axis=1)
    
    if save:
        key = f'{sim}_{str(int(haloid))}'
        filepath = '/home/lonzaric/astro_research/ahdata/ejected_particles.hdf5'
        print(f'Saving {key} ejected particle dataset to {filepath}')
        ejected.to_hdf(filepath, key=key)
        
        filepath = '/home/lonzaric/astro_research/ahdata/cooled_particles.hdf5'
        print(f'Saving {key} cooled particle dataset to {filepath}')
        cooled.to_hdf(filepath, key=key)

        filepath = '/home/lonzaric/astro_research/ahdata/expelled_particles.hdf5'
        print(f'Saving {key} expelled particle dataset to {filepath}')
        expelled.to_hdf(filepath, key=key)
                
        filepath = '/home/lonzaric/astro_research/ahdata/accreted_particles.hdf5'
        print(f'Saving {key} accreted particle dataset to {filepath}')
        accreted.to_hdf(filepath, key=key)
        
        
    print(f'Returning (ejected, cooled, expelled, accreted) datasets.')

    return ejected, cooled, expelled, accreted



# ------------------------------------------- #

def calc_heated(gas, simtime):
    '''
    Overview:
        Quick way to find SN-heated gas particles in general gas distributions.
    
    Inputs:
    -> 'sim', identifier for simulation.
    
    Output:
    -> array of gas particles affected by SN heating.
    '''
    #-----------------------------#
    
    return
    


