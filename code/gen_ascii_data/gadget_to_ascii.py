"""
Code to write ascii snapshots from Gadget binary 2 format.

to-do:
======

- Same mass in satellite and host
- Use bound particles of the host to compute the expansion?

"""

import numpy as np
import sys
sys.path.append("../")
import reading_snapshots as rs



def truncate_halo(pos, rcut):
    r_halo = (pos[:,0]**2 + pos[:,1]**2 + pos[:,2]**2)**0.5
    rcut_index = np.where(r_halo<rcut)[0]
    return pos[rcut_index]

def sample_halo(pos, mass, n_halo_part, npart_sample):
    

    N_random = np.random.randint(0, len(pos), npart_sample)
    	
    mass_fraction = n_halo_part/npart_sample		
    part_mass = mass*mass_fraction
    print('Particle mass factor', mass_fraction)
    print('New particle mass', part_mass)
    return pos[N_random], part_mass

def npart_satellite(pos_sat, pmass_sat, pmass_host):
    """
    Sample satellite galaxies to have the same mass of the host satellite.
    """

    # number of particles in satellite
    init_sat_part = len(pos_sat)
    # Satellite total mass
    sat_tot_mass = pmass_sat*init_sat_part
    # new number of particles of satellite
    n_part_sat = int(sat_tot_mass/pmass_host)
    # sampling satellite particles
    n_part_sample = np.random.randint(0, len(pos_sat), n_part_sat)
    # new particles mass
    new_part_mass = sat_tot_mass/n_part_sat
    return pos_sat[n_part_sample], new_part_mass



def write_snap_txt(path, snap_name, pos, mass):
    np.savetxt(path+snap_name+'.txt', np.array([pos[:,0], pos[:,1], pos[:,2], mass*np.ones(len(pos))]).T)
    return 0

def write_log(halo, sat):
    """
    Printing summary
    """

    print('*********************************************')
    print('Summary:')
    print('Initial number of halo particles: ', halo[0])
    print('Initial halo particle mass: ', halo[1])
    print('Final number of halo particles', halo[2])
    print('Final halo particle mass', halo[3])
    print('Initial number of satellites particles: ', sat[0])
    print('Initial satellites particle mass: ', sat[1])
    print('Final number of satellites particles', sat[2])
    print('Final satellites particle mass', sat[3])
    print('*********************************************')

    return 0


if __name__ == "__main__":
    path = '/media/ngaravito/4fb4fd3d-1665-4892-a18d-bdbb1185a07b1/mwlmc_raw/'
    out_path_MW = '/media/ngaravito/4fb4fd3d-1665-4892-a18d-bdbb1185a07b1/mwlmc_ascii/MW/'
    out_path_LMC = '/media/ngaravito/4fb4fd3d-1665-4892-a18d-bdbb1185a07b1/mwlmc_ascii/LMC/'
    
    #snap_name = 'MWLMC6_2_100M_new_b1_112' 
    #snap_name = 'MWLMC5_100M_new_b1_110' 
    #snap_name = 'MWLMC5_100M_new_b1_110' 
    #snap_name = 'MWLMC5_100M_new_b1_110' 
    #snap_name = "MWLMC3_100M_new_b1_090"
    snap_names = ["MWLMC3_100M_new_b0_090", "MWLMC3_100M_new_b1_091",  
                  "MWLMC4_100M_new_b0_114", "MWLMC4_100M_new_b1_115",  
                  "MWLMC5_100M_new_b1_110", "MWLMC5_100M_new_b0_109",  
                  "MWLMC6_100M_new_b0_2_113", "MWLMC6_100M_new_b1_2_114"]



    n_halo_part = 100000000
    n_part_sample = 100000000
    #n_part_sample_sat = 1000000
    rcut_halo = 400
    sample = 0

    for i in range(len(snap_names)):
        halo = rs.read_snap_coordinates(path, snap_names[i], n_halo_part, com_frame='MW', galaxy='MW')
        # read_snap_coordinates returns pos, vel, pot, mass
        pos_halo_tr = truncate_halo(halo[0], rcut_halo)

        satellite = rs.read_snap_coordinates(path, snap_names[i], n_halo_part, com_frame='sat', galaxy='sat')

        if sample == 1:
            pos_sample, mass_sample = sample_halo(pos_halo_tr, halo[3][0], n_halo_part, n_part_sample)
        elif sample == 0:
            mass_sample = halo[3][0] 
            pos_sample = pos_halo_tr
        pos_sat_em, mass_sat_em =  npart_satellite(satellite[0], satellite[3][0], mass_sample)	
   
        # Outs: 
        out_snap_host = 'MW_{}_{}'.format(int(len(pos_sample)/1E6), snap_names[i])
        out_snap_sat= 'LMC_{}_{}'.format(int(len(pos_sat_em)/1E6), snap_names[i])

        write_log([n_halo_part, halo[3][0], len(pos_sample), mass_sample], [len(satellite[0]), satellite[3][0], len(pos_sat_em), mass_sat_em])
        write_snap_txt(out_path_MW, out_snap_host, pos_sample, mass_sample)
        write_snap_txt(out_path_LMC, out_snap_sat, pos_sat_em, mass_sat_em)
