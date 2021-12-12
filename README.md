# Stellar Feedback-Driven Outflow Analysis
---
#### Author: Leo Lonzarich

#### Director: Professor Charlotte Christensen

#### Institution: Grinnell College, IA

#### Date: 28 Aug. 2021 --> Present 

---

This repository contains the bulk of code used in the analysis of stellar feedback-driven outflows started Fall of 2021. (Currently active). The work here is an extension of computational analysis conducted by Professor Christensen, Hollis Akins, et al., which has been directed to understand galactic star formation quenching in the low-mass, dwarf/satellite regime. This extension of their analysis uses the `DC Justice League Suite` of cosmological, hydrodynamic simulations constructed with Charm N-body Gravity solver (ChaNGa) and prepared by Ferah Munshi (University of Oklahoma) et al. Equipped with these simulations, the computations here take specific aim at gas particles stripped from satellites that have been heated by supernovae (stellar feedback), and is motivated by a desire to understand how such feedback effects, if at all, contribute to satellite quenching, if they aid generally excepted modes of quenching (i.e., ram pressure stripping, tidal stripping), and how they impact satellite morphology.

Note that the datasets used for this study are not included here for sake of their larger sizes, but can be found on Quirm (@grinnell.mathlan) via `/home/lonzaric/astro_research/Stellar_Feedback_Code/SNeData/`.  


<br />


*Disclaimer*. Sections of code utilized here are adapted from Hollis Akins' Justice League Code repository. In the files where such adaptations occur, attributions, along with a permalink to the source file, are given in the header. The permalink to Akins' source repository is given below: 

https://github.com/hollisakins/Justice_League_Code.git

## Satellites

All gas particles treated in this research are sourced from the 19 satellite galaxies (specified by their assigned identifier in each respective simulation) listed as follows: 

- Sim **h148**: 
    - 13, 28, 37, 45, 68, 80, 283, 278, 329

- Sim **h229**:
    - 20, 22, 23, 27, 55
 
- Sim **h242**:
    - 24, 41, 80
    
- Sim **h329**:
    - 33, 137


## Datasets

> `discharged`

- Gas particles that have been removed from the disk of their respective satellite galaxy. This includes both gas moved to halo, and that moved beyond the virial radius. (Properties given for **timestep after discharge event**.)

> `predischarged`

- A collection of the same particles in *discharged*, but instead giving the properties of each particle prior to their discharge event (thereby allowing pre- and post-discharge comparisons).

> `heated`
        
- Gas particles in *discharged* that experienced supernova(e) heating (functionally, their 'cooling' was turned off) in the timestep prior to being discharged. 
- (Note: this dataset can be also be obtained by selecting particles from *discharged* with `sneHeated==True`.)

> `preheated`

- Similar to *predischarged*; Consists of the same particles in *heated*, but with properties of each particle for the timestep prior to discharge.

> `accreted`
        
- Gas particles in the halo or beyond the virial radius of a satellite that are accreted onto the satellite's disk. (Properties given for **timestep following accretion event**.) 
- Note that this includes particles that were previously discharged from the satellite's disk, and those that were not. 

> `reaccreted`

- A subset of *accreted*; Gas particles in the halo or beyond the virial radius of a satellite as a result of a prior discharge event that are reaccreted onto the satellite's disk. (Properties given for **timestep following reaccretion event**.)


## Data Collection Scripts

> `base.py`

- A script containing basic functions used for ubiquitous adminstrative tasks such as obtaining the timesteps for a specified simulation or calling halo keys `{sim}_{haloid}` (e.g. 'h148_13') specifying a satellite galaxy and the simulation it belongs to. Additionally, this script contains a `rootPath` identifier that you should specify as the root directory containing simulation data when running an analysis code in this repository. (Several of the proceeding scripts and notebooks use `base.py`, and will call `base.py` for `rootPath` to streamline directory switching.)

> `particletracking.py`

- An adaption of code from Akins H. by the same name , this script tracks gas particles for a specified simulation and specified satellite. Tracking begins for the first snapshot where the satellite comes within $2R_{vir}$ of its host galaxy, and ends at redshift $z=0$. Any gas particles within the satellite (while it is within this radial limit) for at least one snapshotthis are included.
- Uses as input data: simulation snapshots, `../SNeData/filepaths_haloids.pickle` (Akins H.) and generates a pickle file `../SNeData/iords/{sim}_{haloid}.pickle` and appends to an hdf5 `../SNeData/tracked_particles_v2.hdf5` ('v2' to differentiate from Akins' particle tracking).
- Tracking can be run from the command line via `python particletracking.py h148 13`, where `h148` specifies the simulation and `13` specifies a corresponding satellite. (The bash script `runall.sh` can be used to run multiple instances of particle tracking.) 

> `analysis.py`

- Contains analytical functions to calculated such things as exit angles, used primarily by the compiler when producing datasets for this study.

> `compiler.py`

- This script contains the functions needed to both read and write the .hdf5 datasets used in this treatement. 
- Writing: .hdf5 files are generated in `../SNeData/` from selections of gas particles in `../SNeData/tracked.particles_v2.hdf5`. Useful keys (e.g., `sat_Mvir`, `temp`) are computed and appended here for later use in the analysis notebooks (implemented in the `write_{dataset}.py` scripts). 
- Reading: `read_{dataset}` functions are included, and can be called from, here to read off the .hdf5 files into dataframes.

> `write_discharged.py`, `write_heated.py`, `write_reaccreted.py`

- Each of these scripts runs a `calc_{dataset}` function that is run iteratively over all selected satellites. 
- When run in the terminal (e.g., `python write_discharged.py`), an .hdf5 slice of `../SNeData/tracked_particles_v2.hdf5` with the desired selections of particles. e.g., running `write_discharged.py` has as outputs `../SNeData/predischarged_particles.hdf5`, `../SNeData/discharged_particles.hdf5`, `../SNeData/accreted_particles.hdf5`.


## iPython Notebooks

> `StellarFeedback_p1`

- This notebook includes basic analysis of kinematic et al. properties like exit angles and temperatures for discharged and accreted gas particles. 
- Included here are comparisons of the properties of these particles pre- and post- discharge/accretion, and comparisons of those that were SNe-heated prior to discharged against those that were not.

> `StellarFeedback_p2`

- A follow-up to the prior notebook that includes more rigorous analysis (e.g., using fractional values like $r/R_{vir}$, $v/V_{esc}$ to generate more meaninful results) between gas particles that were and were not SNe-heated just prior to discharge. Additionally, this introduces *reaccreted* particles into the study.

> `StellarFeedback_p3`

- This notebook takes specific aim at studying reaccretion (or lack thereof) of previously discharged particles back onto their respective satellite's disk. This includes looking at reaccretion times and how SNe-heating affects various aspects of reaccretion.
- Analysis here is limited due to time constraints, but will be advanced from Spring 2022.


## Directories

> `../SNeData/`

- Contains all of the data (except for original simulations) used in particle tracking, compiling, and analysis with. While this data is not included here for sake of its size, it can be found in the math.grinnell.edu directory `~/home/lonzaric/astro_research/Stellar_Feedback_Code/`.
- Simulations can be found on Quirm at request.

> `DataAcquisition`

- Contains all `write_{dataset}.py` scripts needed for compiling the specific datasets like *discharged* and *reaccreted* used in the analysis here.


> `plots`

- Contains all plots generated for the final research paper.

> `TestScripts`

- Scripts and notebooks used for troubleshooting kept as reference material.

> `archive`

- Scrapped versions of scripts and notebooks.


<br />


<br />


<br />


*Last rev. 12 Dec. 2021*