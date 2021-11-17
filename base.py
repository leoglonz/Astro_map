# import basic packages
import pynbody
import pandas as pd
import numpy as np
import pickle
import matplotlib.pyplot as plt
import matplotlib as mpl

# define some constants, which should be accessible by any code that imports base.py or analysis.py
hubble =  0.6776942783267969 # hubble constant
age = 13.800797497330507 # age of universe at z=0

# set up matplotlib preferences
mpl.rc('font',**{'family':'serif','monospace':['Palatino']})
mpl.rc('text', usetex=True)
mpl.rcParams.update({'figure.dpi': 200,
                     'font.size': 9,
                     'xtick.direction': 'in',
                     'ytick.direction': 'in',
                     'legend.frameon': False,
                     'figure.constrained_layout.use': True,
                     'xtick.top': True,
                     'ytick.right': True})


# function to plot a median line over some data
def plot_median(ax,x,y,logx=False,logy=False,bins=False,std=False,**kwargs):
    from scipy.stats import binned_statistic
    if logx:
        x = np.log10(x)
        
    if logy:
        y = np.log10(y)
        
    
    condx = ~np.isnan(x) & ~np.isinf(x)
    condy = ~np.isnan(y) & ~np.isinf(y)
    cond = condx & condy
    x, y = x[cond], y[cond]
        
    if type(bins)==bool:
        bins = np.linspace(np.min(x), np.max(x), 10)
    if type(bins)==int:
        bins = np.linspace(np.min(x), np.max(x), bins)
    
    # calculate median
    median, bins, binnumber = binned_statistic(x,y,bins=bins,statistic='median')
    bc = 0.5*(bins[1:]+bins[:-1])
    
    if logx:
        bc = np.power(10,bc)
        
    if std:
        std, bins, binnumber = binned_statistic(x,y,bins=bins,statistic='std')
        if 'color' in kwargs:
            mycolor = kwargs.get('color')
        else:
            mycolor = 'tab:blue'
            
        ymin, ymax = median-std, median+std
        if logy:
            ymin, ymax = np.power(10,ymin), np.power(10,ymax)
            
        ax.fill_between(bc,ymin, ymax, fc=mycolor, ec=None, alpha=0.15)

        
    if logy:
        median = np.power(10,median) 
        
    ax.plot(bc, median, **kwargs)
    
setattr(mpl.axes.Axes, "plot_median", plot_median)



# define functions for basic data manipulation, importing, etc. used by everything
def get_stored_filepaths_haloids(sim,z0haloid):
    # get snapshot paths and haloids from stored file
    with open('/home/lonzaric/astro_research/Stellar_Feedback_Code/SNeData/filepaths_haloids.pickle','rb') as f:
        d = pickle.load(f)
    try:
        filepaths = d['filepaths'][sim]
    except KeyError:
        print("sim must be one of 'h148','h229','h242','h329'")
        raise
    try:
        haloids = d['haloids'][sim][z0haloid]
        h1ids = d['haloids'][sim][1]
    except KeyError:
        print('z0haloid not found, perhaps this is a halo that has no stars at z=0, and therefore isn\'t tracked')
        raise
    return filepaths,haloids,h1ids
    
# timesteps data
def read_timesteps(sim):
    '''Function to read in the data file which contains quenching and infall times'''
    data = []
    with open(f'/home/lonzaric/astro_research/Stellar_Feedback_Code/SNeData/timesteps_data/{sim}.data', 'rb') as f:
        while True:
            try:
                data.append(pickle.load(f,encoding='latin1'))
            except EOFError:
                break

    data = pd.DataFrame(data)
    return data

# timescales (quenching timescales, derived from timesteps)
def read_timescales():
    '''Function to read in the data file which contains quenching and infall times'''
    data = []
    with open('/home/lonzaric/astro_research/Stellar_Feedback_Code/SNeData/QuenchingTimescales.data', 'rb') as f:
        while True:
            try:
                data.append(pickle.load(f,encoding='latin1'))
            except EOFError:
                break

    data = pd.DataFrame(data)
    return data

# infall properties (properties of satellites at t_infall, used in Figure 6 of Akins et al. 2021)
def read_infall_properties():
    '''Function to read in the data file with quenching timescales and satellite properties at infall.'''
    data = []
    with open(f'/home/lonzaric/astro_research/Stellar_Feedback_Code/SNeData/QuenchingTimescales_InfallProperties.data','rb') as f:
        while True:
            try: 
                data.append(pickle.load(f))
            except EOFError:
                break
            
    data = pd.DataFrame(data)
    data['timescale'] = data.tinfall - data.tquench
    
    return data

# determines the snapshot at which to start tracking (first snapshot where satellite is within 2 Rvir of host)
def get_snap_start(sim,z0haloid):
    print(f'\t {sim}-{z0haloid}: Getting starting snapshot (dist = 2 Rvir)')
    filepaths,haloids,h1ids = get_stored_filepaths_haloids(sim,z0haloid)
    ts = read_timesteps(sim)
    ts = ts[ts.z0haloid == z0haloid]

    dist = np.array(ts.h1dist, dtype=float)
    time = np.array(ts.time, dtype=float)
    ti = np.min(time[dist <= 2])

    for i,f in enumerate(filepaths):
        s = pynbody.load(f)
        t = float(s.properties['time'].in_units('Gyr'))
        if t<ti:
            snap_start = i
            break
        else: 
            continue
    print(f'\t {sim}-{z0haloid}: Start on snapshot {snap_start}, {filepaths[snap_start][-4:]}') # go down from there!
    return snap_start