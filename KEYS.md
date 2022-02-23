## The Comprehensive, Non-exhaustive List of Incorporated Properties
---

Concerning the gas particles processed by Hollis Akins et al. in the near-mint `200_bkgdens Justice League Simulations`: explanations for assorted properties (`keys`) and their intrinsic units denoted in a parenthetical endnote for each.
<br><br>
There are three primary categories of properties attached to particle data:
- Properties relative to satellite;
- Properties relative to host;
- Properties of a particle's satellite galaxy relative to the host galaxy.

Where used, the center of a satellite/host galaxy is defined by its center of mass (COM).

---



## The Quick Version: 
### In General (unless where otherwise noted)...
- *Time* in Gigayears $(Gyrs)$.
- *Distances* in kiloparsecs $(kpc)$.
- *Velocities* in kilometers per second $(km \ s^{-1})$.
- *Masses* in solar masses $(M_{sol})$.
- *Gravitational Constant* G in ($kpc\ km^2\ M_{sol}^{-1}\ s^{-2})$



## The Long Version:
### Satellite-oriented Properties: (satellite center == COM)

- `r`: radial distance from satellite. (kpc)
- `r_per_Rvir`:  `r` per viral radius of satellite (sim. viral radius scaled by $$a/Hubble_constant$$).
- `x`, `y`, `z`: coordinates relative to satellite. (kpc)
- `satRvir`: virial radius of satellite. (kpc)
- `a`: Justice League sim. constant.
- `vx`, `vy`, `vz`: component velocities rel. satellite. (km s-1)
- `v`:  speed rel. to satellite. (km s-1)


### Host-oriented Properties: (host centre == COM)

- `r_rel_host`: radial position from host. (kpc)
- `r_rel_host_per_Rvir`: `r_rel_host` per virial radius of host (sim. viral radius scaled by $$1/a * Hubble constant$$).
- `x_rel_host`, `y_rel_host`, `z_rel_host`: coord. positions rel. to host. (kpc)
- `hostRvir`: virial radius of host. (kpc)
- `vx_rel_host`,`vy_rel_host`, `vz_rel_host`: velocities rel. to host. (km s-1)
- `v_rel_host`: speed rel. host. (km s-1)


### Properties of Satellite Relative to Host: (sat./host centers == COM)

- `sat_Xc`,`sat_Yc`, `satZc`: coords of sat. relative to center of simulation snapshot (scaled by $1/a * Hubble constant$). (kpc)
- `sat_vx`,`sat_vy`, `sat_vz`: velocities of sat. relative to cener of sim. (km s-1)
- `host_Xc`, `host_Yc`, `host_Zc`: coords of host relative to center of sim. (kpc)
- `host_vx`, `host_vy`, `host_vz`: velocities of host relative to center of sim. (km s-1) 




*Last rev. 7 Feb. 2022*