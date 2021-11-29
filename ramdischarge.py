# Contains all functions required to calculate and compute ram pressures for the gas particle datasets.
#
# _______________________________________________________________________________________________
# Ram pressure calculations credit to Hollis Akins 2021;
# Github permalink: https://github.com/hollisakins/Justice_League_Code/blob/
#                    e049137edcfdc9838ebb3cf0fcaa4ee46e977cec/Analysis/RamPressure/rampressure.py
# _______________________________________________________________________________________________
# Last revised: 29 Nov. 2021

import sys
import tqdm
import os
# import fsps

from base import *
from compiler import *
from analysis import *



def calc_ram_pressure(sim, z0haloid, filepaths, haloids, h1ids):  
    '''
    -> Ram pressure calculations for desired sets of gas particles..
    '''
    #--------------------------------#
    output_tot = pd.DataFrame()
    
    #print('Starting calculations...')
    for f,haloid,h1id in tqdm.tqdm(zip(filepaths,haloids,h1ids),total=len(filepaths)):
        # load simulation
        s = pynbody.load(f)
        s.physical_units()
        h = s.halos()
        sat = h[haloid]
        host = h[h1id]
        
        # save time t and scale factor a
        t = float(s.properties['time'].in_units('Gyr'))
        a = float(s.properties['a'])

        output = pd.DataFrame()
        output['t'] = [t]
        output['a'] = [a]
        
        # positions and velocities
        r_sat = np.array([sat.properties[k]/hubble*a for k in ['Xc','Yc','Zc']])
        r_host = np.array([host.properties[k]/hubble*a for k in ['Xc','Yc','Zc']])
        r_rel = r_sat - r_host
        h1dist = np.linalg.norm(r_rel)
        output['h1dist'] = [h1dist]
        print(f'\n\t {sim}-{z0haloid}: Distance from host = {h1dist:.2f} kpc')
        
        v_sat = np.array([sat.properties[k] for k in ['VXc','VYc','VZc']])
        v_host = np.array([host.properties[k] for k in ['VXc','VYc','VZc']])
        v_rel = v_sat - v_host
        v_rel_mag = np.linalg.norm(v_rel)
        print(f'\t {sim}-{z0haloid}: Relative velocity = {v_rel_mag:.2f} km/s')

        # nearest neighbor distance (topic of ongoing investigation)
        nDM = len(sat.dm)
        Xcs, Ycs, Zcs = np.array([]), np.array([]), np.array([])
        for halo in h:
            if len(halo.dm) > nDM*0.1:
                r = np.array([halo.properties[k]/hubble*a for k in ['Xc','Yc','Zc']])
                if not (r==r_sat).all():
                    Xcs = np.append(Xcs, r[0])
                    Ycs = np.append(Ycs, r[1])
                    Zcs = np.append(Zcs, r[2])
        
        
        r_others = np.array([Xcs, Ycs, Zcs]).T
        dists = np.linalg.norm(r_others - r_sat, axis=1)
        print(f'\t {sim}-{z0haloid}: dNN = {np.min(dists):.1f} kpc')
        output['dNN'] = [np.min(dists)]
        
        # basic galaxy properties
        M_star = np.sum(sat.s['mass'].in_units('Msol'))
        M_gas = np.sum(sat.g['mass'].in_units('Msol'))
        rvir = sat.properties['Rvir']/hubble*a
        h1rvir = host.properties['Rvir']/hubble*a

        output['M_star'] = [M_star]
        output['M_gas'] = [M_gas]
        output['satRvir'] = [rvir]
        output['hostRvir'] = [h1rvir]
        print(f'\t {sim}-{z0haloid}: Satellite M_gas = {M_gas:.1e} Msun')

        # simple ram pressure calculations: calculate rho_CGM from spherical density profile
        pynbody.analysis.halo.center(host)
        pg = pynbody.analysis.profile.Profile(s.g, min=0.01, max=2*h1dist, ndim=3)
        rbins = pg['rbins']
        density = pg['density']

        rho_CGM = density[np.argmin(np.abs(rbins-h1dist))]
        Pram = rho_CGM * v_rel_mag * v_rel_mag
        output['vel_CGM'] = [v_rel_mag]
        output['rho_CGM'] = [rho_CGM]
        output['Pram'] = [Pram]
        print(f'\t {sim}-{z0haloid}: Simple v_rel = {v_rel_mag:.1f}')
        print(f'\t {sim}-{z0haloid}: Simple rho_CGM = {rho_CGM:.1e}')
        print(f'\t {sim}-{z0haloid}: Simple P_ram = {Pram:.1e}')


        # advanced ram pressure calculations: calculate rho, vel from cylinder in front of satellite
        
        # code below is adapted from pynbody.analysis.angmom.sideon()
        # transform the snapshot so that the vector 'vel' points in the +y direction
        top = s
        print(f'\t {sim}-{z0haloid}: Centering positions')
        cen = pynbody.analysis.halo.center(sat, retcen=True)
        tx = pynbody.transformation.inverse_translate(top, cen)
        print(f'\t {sim}-{z0haloid}: Centering velocities')
        vcen = pynbody.analysis.halo.vel_center(sat, retcen=True)
        tx = pynbody.transformation.inverse_v_translate(tx, vcen)
        
        # try to get vel from gas particles, but if there are no gas particles, use stars
        print(f'\t {sim}-{z0haloid}: Getting velocity vector')
        try:
            vel = np.average(sat.g['vel'], axis=0, weights=sat.g['mass'])
        except ZeroDivisionError:
            vel = np.average(sat.s['vel'], axis=0, weights=sat.s['mass'])
            
        vel_host = np.average(host.g['vel'], axis=0, weights=host.g['mass'])
        vel -= vel_host

        print(f'\t {sim}-{z0haloid}: Transforming snapshot')
        trans = vec_to_xform(vel)
        tx = pynbody.transformation.transform(tx, trans)
        
        # define cylinder size and filter out those particles
        radius = 0.5*rvir
        height = 0.75 * radius
        center = (0, rvir + height/2, 0)
        wind_filt = pynbody.filt.Disc(radius, height, cen=center)
        env = s[wind_filt].g
        print(f'\t {sim}-{z0haloid}: Identified {len(env)} gas particles to calculate wind properties')
        output['n_CGM'] = [len(env)]

        # try to calculate CGM properties, but if you can't then set rho, vel to 0 (i.e. no gas particles)
        try:
            vel_CGM = np.linalg.norm(np.average(env['vel'],axis=0,weights=env['mass'])) # should be in units of Msun kpc**-3
            rho_CGM = np.average(env['rho'], weights=env['mass']) # should be in units of
            std_rho_CGM = np.std(env['rho'])
            std_vel_CGM = np.std(np.linalg.norm(env['vel'], axis=1))
        except ZeroDivisionError:
            vel_CGM, rho_CGM = 0, 0
            std_vel_CGM, std_rho_CGM = 0, 0
            
        Pram = rho_CGM * vel_CGM * vel_CGM # overall units should be Msun kpc**-3 km s**-1
        
        output['vel_CGM_adv'] = [vel_CGM]
        output['rho_CGM_adv'] = [rho_CGM]
        output['std_vel_CGM'] = [std_vel_CGM]
        output['std_rho_CGM'] = [std_rho_CGM]
        output['Pram_adv'] = [Pram]

        print(f'\t {sim}-{z0haloid}: Advanced vel_CGM = {vel_CGM:.2f}')
        print(f'\t {sim}-{z0haloid}: Advanced rho_CGM = {rho_CGM:.1e}')
        print(f'\t {sim}-{z0haloid}: Advanced P_ram = {Pram:.1e}')

        # restoring pressure calculations
        # try to center the satellite. if you can't, then that means Prest = 0 (i.e. Mgas=0)
        try:
            pynbody.analysis.halo.center(sat)
            calc_rest = True
        except:
            calc_rest = False
            Prest = 0.
            SigmaGas = 0.
            dphidz = 0.
        
        if calc_rest:
            p = pynbody.analysis.profile.Profile(s.g, min=0.01, max=rvir, ndim=3)
            percent_enc = p['mass_enc']/M_gas
            rhalf = np.min(p['rbins'][percent_enc > 0.5])
            SigmaGas = M_gas / (2*np.pi*rhalf**2)
            Rmax = sat.properties['Rmax']/hubble*a
            Vmax = sat.properties['Vmax']
            dphidz = Vmax**2 / Rmax
            Prest = dphidz * SigmaGas
        
        print(f'\t {sim}-{z0haloid}:  Prest = {Prest:.1e}')

        output['Prest'] = [Prest]
        output['SigmaGas'] = [SigmaGas]
        output['dphidz'] = [dphidz]
    
        
        # sfr calculations: use FSPS to calculate SFRs from formation masses of stars, not current masses
        star_masses = np.array(sat.s['mass'].in_units('Msol'),dtype=float)
        star_metals = np.array(sat.s['metals'], dtype=float)
        star_ages = np.array(sat.s['age'].in_units('Myr'),dtype=float)
        size = len(star_ages)
        
        # construct simple stellar population with fsps
        fsps_ssp = fsps.StellarPopulation(sfh=0,zcontinuous=1, imf_type=2, zred=0.0, add_dust_emission=False)
        solar_Z = 0.0196
        
        star_masses = star_masses[star_ages <= 100]
        star_metals = star_metals[star_ages <= 100]
        star_ages = star_ages[star_ages <= 100]
        print(f'\t {sim}-{z0haloid}: performing FSPS calculations on {len(star_masses)} star particles (subset of {size} stars)')
        
        if len(star_masses)==0:
            SFR = 0
        else:
            massform = np.array([])
            for age, metallicity, mass in zip(star_ages, star_metals, star_masses):
                fsps_ssp.params['logzsol'] = np.log10(metallicity/solar_Z)
                mass_remaining = fsps_ssp.stellar_mass
                massform = np.append(massform, mass / np.interp(np.log10(age*1e9), fsps_ssp.ssp_ages, mass_remaining))

            SFR = np.sum(massform)/100e6
            
        output['SFR'] = [SFR]
        output['sSFR'] = [SFR/M_star]
        print(f'\t {sim}-{z0haloid}: sSFR = {SFR/M_star:.2e} yr**-1')
        
        output_tot = pd.concat([output_tot, output])

    return output_tot



def read_ram_pressure(sim, haloid):
    '''
    Function to read in the ram pressure dataset, merge it with particle and flow information, and return a dataset containing 
    rates of gas flow in addition to ram pressure information.
    '''
    #--------------------------------#

    # loading ram pressure data for specified simulation and haloid.
    path = f'{rootPath}Stellar_Feedback_Code/SNeData/ram_pressure.hdf5'
    key = f'{sim}_{haloid}'
    data = pd.read_hdf(path, key=key)
    
    # converting data to numpy arrays (i.e. remove pynbody unit information) and calculating ratio
    data['Pram_adv'] = np.array(data.Pram_adv,dtype=float)
    data['Pram'] = np.array(data.Pram,dtype=float)
    data['Prest'] = np.array(data.Prest,dtype=float)
    data['ratio'] = data.Pram_adv / data.Prest
    dt = np.array(data.t)[1:] - np.array(data.t)[:-1]
    dt = np.append(dt[0],dt)
    data['dt'] = dt
    
    # loading timescale information to add quenching time and quenching timescale (tau).
    timescales = read_timescales()
    ts = timescales[(timescales.sim==sim)&(timescales.haloid==haloid)]
    data['tau'] = ts.tinfall.iloc[0] - ts.tquench.iloc[0]    
    data['tquench'] = age - ts.tquench.iloc[0]   

    # loading discharged particle data.
    predischarged, discharged, accreted, preheated, heated = read_all_discharged();

    # Mgas_div is the gas mass we divide by when plotting rates. this is the gas mass 1 snapshot past.
    Mgas_div = np.array(data.M_gas,dtype=float)
    Mgas_div = np.append(Mgas_div[0], Mgas_div[:-1])
    data['Mgas_div'] = Mgas_div
    
    # load in particle data.
    particles = read_tracked_particles(sim,haloid)
    # m_disk = 0 if particle is not in disk, = particle mass if it is. This allows us to compute total mass in the disk.
    particles['m_disk'] = np.array(particles.mass,dtype=float)*np.array(particles.sat_disk,dtype=int)
    particles['m_SNeaff'] = np.array(particles.mass,dtype=float)*np.array(particles.coolontime > particles.time, dtype=int)
    
    # group particles data by unique times and sum the mass of particles that are SNe affected, to get total mass.
    data = pd.merge_asof(data, particles.groupby(['time']).m_SNeaff.sum().reset_index(), left_on='t', right_on='time')
    data = data.rename(columns={'m_SNeaff':'M_SNeaff'})
    
    # group particle data by unique times and sum the mass of particles that are in the disk, to get total mass.
    data = pd.merge_asof(data, particles.groupby(['time']).m_disk.sum().reset_index(), left_on='t', right_on='time')
    data = data.rename(columns={'m_disk':'M_disk'})
    
    # analagous to Mgas_div above.
    Mdisk_div = np.array(data.M_disk,dtype=float)
    Mdisk_div = np.append(Mdisk_div[0], Mdisk_div[:-1])
    data['Mdisk_div'] = Mdisk_div
    
    # 1) fetching rates of predischarged gas.
    data = pd.merge_asof(data, predischarged.groupby(['time']).mass.sum().reset_index(), left_on='t', right_on='time')
    data = data.rename(columns={'mass':'M_predischarged'}) # mass ejected in that snapshot
    data['Mdot_predischarged'] = data.M_predischarged / data.dt # rate of mass ejection 
    data['Mdot_predischarged_by_Mgas'] = data.Mdot_predischarged / Mgas_div # rate of ejection divided by M_gas
    data['Mdot_predischarged_by_Mdisk'] = data.Mdot_predischarged / Mdisk_div # rate of ejection divided by M_disk

    # 2) fetching rates of all discharged gas.
    data = pd.merge_asof(data, discharged.groupby(['time']).mass.sum().reset_index(), left_on='t', right_on='time')
    data = data.rename(columns={'mass':'M_discharged'}) 
    data['Mdot_discharged'] = data.M_discharged / data.dt 
    data['Mdot_discharged_by_Mgas'] = data.Mdot_discharged / Mgas_div 
    data['Mdot_discharged_by_Mdisk'] = data.Mdot_discharged / Mdisk_div 

    # 3) fetching for gas reaccreted onto satellite disk.
    data = pd.merge_asof(data, accreted.groupby(['time']).mass.sum().reset_index(), left_on='t', right_on='time')
    data = data.rename(columns={'mass':'M_accreted'})
    data['Mdot_accreted'] = data.M_accreted / data.dt
    data['Mdot_accreted_by_Mgas'] = data.Mdot_accreted / Mgas_div
    data['Mdot_accreted_by_Mdisk'] = data.Mdot_accreted / Mdisk_div
    
    # 4) fetching for predischarged heated gas.
    data = pd.merge_asof(data, preheated.groupby(['time']).mass.sum().reset_index(), left_on='t', right_on='time')
    data = data.rename(columns={'mass':'M_preheated'})
    data['Mdot_preheated'] = data.M_preheated / data.dt
    data['Mdot_preheated_by_Mgas'] = data.Mdot_preheated / Mgas_div
    data['Mdot_preheated_by_Mdisk'] = data.Mdot_preheated / Mdisk_div
    
    # 5) fetching for discharged heated gas.
    data = pd.merge_asof(data, heated.groupby(['time']).mass.sum().reset_index(), left_on='t', right_on='time')
    data = data.rename(columns={'mass':'M_heated'})
    data['Mdot_heated'] = data.M_accreted / data.dt
    data['Mdot_heated_by_Mgas'] = data.Mdot_heated / Mgas_div
    data['Mdot_heated_by_Mdisk'] = data.Mdot_heated / Mdisk_div

    # overall rate of gas-loss
    dM_gas = np.array(data.M_gas,dtype=float)[1:] - np.array(data.M_gas,dtype=float)[:-1]
    dM_gas = np.append([np.nan],dM_gas)
    data['Mdot_gas'] = dM_gas / np.array(data.dt)
    
    # rate of gas-loss from the disk
    dM_disk = np.array(data.M_disk,dtype=float)[1:] - np.array(data.M_disk,dtype=float)[:-1]
    dM_disk = np.append([np.nan],dM_disk)
    data['Mdot_disk'] = dM_disk / np.array(data.dt)
    
    data['key'] = key
    
    # fraction of the inital gas mass still remaining in the satellite
    M_gas_init = np.array(data.M_gas)[np.argmin(data.t)]
    data['f_gas'] = np.array(data.M_gas)/M_gas_init
    
    return data
 
    
    
def read_all_ram_pressure():
    '''
    -> Returns workable dataframes containing ram pressures for specified set of gas particles.
    '''
    #--------------------------------#    
    
    data_all = pd.DataFrame();
    
    keys = ['h148_13','h148_28','h148_37','h148_45','h148_68','h148_80','h148_283',
            'h148_278','h148_329','h229_20','h229_22','h229_23','h229_27','h229_55',
            'h242_24','h242_41','h242_80','h329_33','h329_137']
    
    i = 1
    for key in keys:
        print(i, end=' ')
        i += 1
        sim = key[:4]
        haloid = int(key[5:])
        data = read_ram_pressure(sim, haloid)
        data_all = pd.concat([data_all,data])  
    
    return data_all



if __name__ == '__main__':
    sim = str(sys.argv[1])
    z0haloid = int(sys.argv[2])
    
    snap_start = get_snap_start(sim,z0haloid)
    filepaths, haloids, h1ids = get_stored_filepaths_haloids(sim,z0haloid)

    # fix the case where the satellite doesn't have merger info prior to
    if len(haloids) < snap_start:
        snap_start = len(haloids)
        raise Exception('Careful! You may have an error since the satellite doesnt have mergertree info out to the time where you want to start. This case is untested')
    
    if len(haloids) >= snap_start:
        filepaths = np.flip(filepaths[:snap_start+1])
        haloids = np.flip(haloids[:snap_start+1])
        h1ids = np.flip(h1ids[:snap_start+1])

    output_tot = calc_ram_pressure(sim, z0haloid, filepaths, haloids, h1ids)
    output_tot.to_hdf('../../Data/ram_pressure_dsrg.hdf5',key=f'{sim}_{z0haloid}')
