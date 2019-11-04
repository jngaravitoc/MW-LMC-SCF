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
    header = 'x[kpc], y[kpc], z[kpc], vx[km/s], vy[km/s], vz[km/s]'
    np.savetxt(filename, data, header=header)
    
def distance_cut(pos, vel, rmin, rmax):
    r = np.sqrt(pos[:,0]**2+pos[:,1]**2+pos[:,2]**2)
    rcut = np.where((r<rmax) & (r>rmin))
    return pos[rcut], vel[rcut]

if __name__ == "__main__":


    path = '/home/xzk/work/github/MW_anisotropy/code/test_snaps/'
    snapshots = ['MWLMC6_2_100M_new_b1_112','MWLMC6_2_100M_new_b0_111',
            'MWLMC5_100M_new_b1_110', 'MWLMC5_100M_new_b0_109',
            'MWLMC4_100M_new_b1_115', 'MWLMC4_100M_new_b0_114',
            'MWLMC3_100M_new_b1_091', 'MWLMC3_100M_new_b0_090']
    N_halo_part = 100000000
    filenames = ['MWLMC6_b1_dm_halo_020_200.txt','MWLMC6_b0_dm_halo_020_200.txt', 
            'MWLMC5_b1_dm_halo_020_200.txt', 'MWLMC5_b0_dm_halo_020_200.txt',
            'MWLMC4_b1_dm_halo_020_200.txt', 'MWLMC4_b0_dm_halo_020_200.txt',
            'MWLMC3_b1_dm_halo_020_200.txt', 'MWLMC3_b0_dm_halo_020_200.txt']
    for i in range(0, len(filenames)):
        MW_particles = reading_snapshots.read_MW_snap_com_coordinates(path, snapshots[i], LMC=True, N_halo_part=N_halo_part, pot=True)
        pos_MW = MW_particles[0]
        vel_MW = MW_particles[1]
        pos_cut, vel_cut = distance_cut(pos_MW, vel_MW, 20, 200)
        write_snapshot(filenames[i], pos_cut, vel_cut)
