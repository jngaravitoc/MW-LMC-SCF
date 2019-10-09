"""
Code to extract information from the DM halos in the N-body simulations.

Author: github/jngaravitoc
Date: 10/08/19

Dependecies:
============

reading_snapshots.py
"""
import numpy as np
import reading_snapshots


def write_snapshot(filename, pos, vel):
    data = np.array([pos[:,0], pos[:,1], pos[:,2], vel[:,0], vel[:,1], vel[:,2]]).T
    np.savetxt(filename, data)
    
def distance_cut(pos, vel, rmin, rmax):
    r = np.sum(pos**2, axis=-1)
    rcut = np.where((r<rmax) & (r>rmin))
    return pos[rcut], vel[rcut]

if __name__ == "__main__":


    path = '../../MW_anisotropy/code/test_snaps/'
    snaphot = 'MWLMC5_100M_new_b1_110'
    N_halo_part = 100000000

    MW_particles = reading_snapshots.read_MW_snap_com_coordinates(path, snapshot, LMC=True, N_halo_part=N_halo_part, pot=False)
    pos_MW = MW_particles[0]
    vel_MW = MW_particles[1]
    pos_cut, vel_cut = distance_cut(pos_MW, vel_MW, 20, 40)
    print(np.len(pos_cut))
