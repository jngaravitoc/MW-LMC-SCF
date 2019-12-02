"""
TODO:

	import all functions properly from LMC_bounded and gadget_to_ascii
	include all input parameters
	implement all optional outputs
	implement checks!

"""

import numpy as np
import sys
import LMC_bounded
import gadget_to_ascii

if __name__ == "__main__":
    path = "../../../MW_anisotropy/code/test_snaps/"
    out_path_MW = './'
    out_path_LMC = './'

    #snap_names = ["MWLMC3_100M_new_b0_090",
    #			  "MWLMC3_100M_new_b1_091",
    #			  "MWLMC4_100M_new_b0_114",
    #             "MWLMC4_100M_new_b1_115",
    #			  "MWLMC5_100M_new_b1_110",
    #			  "MWLMC5_100M_new_b0_109",
    #			  "MWLMC6_100M_new_b0_2_113",
    #			  "MWLMC6_100M_new_b1_2_114"]

    in_path = sys.argv[1]
    snapname = sys.argv[2]
    out_name = sys.argv[3]
    
    nmax = int(sys.argv[4])
    lmax = int(sys.argv[5])
    rs = float(sys.argv[6])
    #n_halo_part = 100000000
    n_halo_part = int(sys.argv[7])
	
    n_part_sample = 100000000
    #n_part_sample_sat = 1000000
    rcut_halo = 400
    sample = 0
    sample_lmc = 0
    #for i in range(0, len(snap_names)):
    for i in range(init_snap, final_snap)
        halo = rs.read_snap_coordinates(in_path, snapname+"{}".format(i), n_halo_part, com_frame='MW', galaxy='MW')
        # read_snap_coordinates returns pos, vel, pot, mass
        pos_halo_tr, vel_halo_tr, mass_tr, ids_tr = truncate_halo(halo[0], halo[1], halo[3], halo[4], rcut_halo)

        satellite = rs.read_snap_coordinates(in_path, snapname+"{}".format(i), n_halo_part, com_frame='sat', galaxy='sat')

        print("**************************")
        print(snapname+"{}".format(i))
        #print(pos_cm, vel_cm)
        print("**************************")

        pos_sat_tr, vel_sat_tr, mass_sat_tr, ids_sat_tr = truncate_halo(satellite[0], satellite[1], satellite[3], satellite[4], rcut_halo)
        pos_sat_em, vel_sat_em, mass_sat_em, ids_sat_em = npart_satellite(pos_sat_tr, vel_sat_tr, ids_sat_tr, mass_sat_tr[0], mass_tr[0])

        # Outs: 
        out_snap_host = 'MW_{}_{}'.format(int(len(pos_halo_tr)/1E6), snapname+"{}".format(i))
        out_snap_sat= 'LMC_{}_{}'.format(int(len(pos_sat_em)/1E6), snapname+"{}".format(i))

        #write_log([n_halo_part, halo[3][0], len(pos_sample), mass_sample], [len(pos_sat_tr[0]), satellite[3][0], len(pos_sat_em), mass_sat_em])
        write_snap_txt(out_path_MW, out_snap_host, pos_halo_tr, vel_halo_tr, mass_tr, ids_tr)
        write_snap_txt(out_path_LMC, out_snap_sat, pos_sat_em, vel_sat_em, mass_sat_em, ids_sat_em)
        
		#write_snap_txt(out_path_LMC, out_snap_sat, satellite[0], satellite[1], satellite[3], satellite[4])
	
		## Satellite bound particles
		#pos, vel, mass, ids = reading_particles(snapname)
		armadillo = find_bound_particles(pos_sat_em, vel_sat_em, mass_sat_em, 
										 ids_sat_em, rs, nmax, lmax)
		print('Bound particles computed')

		pos_bound = armadillo[0]
		vel_bound = armadillo[1]
		N_bound =  armadillo[2]
		ids_bound =  armadillo[3]
		pos_unbound = armadillo[4]
		vel_unbound = armadillo[5]
		ids_unbound = armadillo[6]

		lmc_bound = np.array([pos_bound[:,0], pos_bound[:,1], pos_bound[:,2],
		                      vel_bound[:,0], vel_bound[:,1], vel_bound[:,2],
		                      ids_bound]).T

		lmc_unbound = np.array([pos_unbound[:,0], pos_unbound[:,1], pos_unbound[:,2],
		                        vel_unbound[:,0], vel_unbound[:,1], vel_unbound[:,2],
		                        ids_unbound]).T
		print('Combining satellite unbound particles with host particles')
		
		pos_host_sat = np.vstack((pos_halo_tr, pos_unbound))		
		vel_host_sat = np.vstack((vel_halo_tr, vel_unbound))
		ids_host_sat = np.hstack((ids_halo_tr, ids_unbound))
		## TODO: include out_path
		np.savetxt(out_name, lmc_bound)
		print('Done writing snapshot with satellite bounded particles')
		np.savetxt("unbound"+out_name, lmc_unbound)
		print('Done writing snapshot with satellite unbounded particles')

