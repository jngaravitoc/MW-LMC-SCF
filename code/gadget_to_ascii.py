"""
Code to write ascii snapshots from Gadget binary 2 format.

to-do:
======

- Use bound particles of the host to compute the expansion?

"""

import numpy as np
import pygadgetreader 
import reading_snapshots


def read_mw_com_snapshot(path, snap_name, n_halo_part):
        pos, vel, pot, ids = reading_snapshots.read_MW_snap_com_coordinates(path, snap_name, LMC=True,
                                                                       N_halo_part=n_halo_part, pot=True)
        mass =  pygadgetreader.readsnap(path+snap_name, 'mass', 'dm')
        all_ids =  pygadgetreader.readsnap(path+snap_name, 'pid', 'dm')
        mw_mass = np.where(all_ids==ids[0])
        return pos, mass[mw_mass]


def read_sat_com_snapshot(path, snap_name, n_halo_part):
        pos, vel, lmc_ids = reading_snapshots.read_satellite_snap_com_coordinates(path+snap_name, LMC=True,
                                                                    N_halo_part=n_halo_part, pot=True)
        mass =  pygadgetreader.readsnap(path+snap_name, 'mass', 'dm')
        ids =  pygadgetreader.readsnap(path+snap_name, 'pid', 'dm')
        lmc_mass = np.where(ids==lmc_ids[0])
        
        return pos, mass[lmc_mass], len(ids)


def truncate_halo(pos, rcut):
	r_halo = (pos[:,0]**2 + pos[:,1]**2 + pos[:,2]**2)**0.5
	rcut_index = np.where(r_halo<rcut)[0]
	return pos[rcut_index]

def sample_halo(pos, mass, n_halo_part, npart_sample):
	assert(npart_sample<n_halo_part)

	N_random = np.random.randint(0, len(pos), npart_sample)
	
	mass_fraction = n_halo_part/npart_sample		
	part_mass = mass*mass_fraction
	print('Particle mass factor', mass_fraction)
	print('New particle mass', part_mass)
	return pos[N_random], mass_fraction


def write_snap_txt(snap_name, pos, mass):
	np.savetxt(snap_name+'.txt', np.array([pos[:,0], pos[:,1], pos[:,2], mass*np.ones(len(pos))]).T)


if __name__ == "__main__":
        path = '../../MW_anisotropy/code/test_snaps/'
        snap_name = 'MWLMC5_100M_new_b1_110' 
        #snap_name = 'MWLMC5_100M_new_b1_110' 
        #snap_name = 'MWLMC5_100M_new_b1_110' 
        #snap_name = 'MWLMC5_100M_new_b1_110' 
        n_halo_part = 100000000
        n_part_sample = 1000000
        n_part_sample_sat = 1000000
        rcut_halo = 400
        pos_halo, mass_halo = read_mw_com_snapshot(path, snap_name, n_halo_part)
        pos_halo_tr = truncate_halo(pos_halo, rcut_halo)
        pos_satellite, mass_satellite , n_sat = read_sat_com_snapshot(path, snap_name, n_halo_part)
        print(mass_halo, mass_satellite)
	
        pos_sample, mass_sample = sample_halo(pos_halo_tr, mass_halo, n_halo_part, n_part_sample)
        pos_sat_sample, mass_sat_sample = sample_halo(pos_satellite, mass_satellite, n_sat, n_part_sample_sat)
	#write_snap_txt(path+snap_name, pos_sample, mass_sample)
	#write_snap_txt(path+snap_name_sat, pos_sat_sample, mass_sat_sample)

	
