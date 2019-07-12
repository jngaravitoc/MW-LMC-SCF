import numpy as np
import biff
import coefficients_smoothing
import multiprocessing
from joblib import Parallel, delayed


def kpq_sn_estimator(i):
    pos = np.loadtxt('./MW_100M_b1_dm_part_1e6_300.txt')

    N_particles = len(pos)
    n_batches = 100
    particles_in_batch = int(N_particles/n_batches)


    sn_range = np.arange(0, 11, 0.2)
    Hp = np.zeros(len(sn_range))
    Kpq = np.zeros(len(sn_range))
    Kpq_rho = np.zeros(len(sn_range))
    N_part = 990000 # len(pos_grid)
    m_p = 1/N_part
    rho_factor = 1e4/N_part
    m_p_sim = 1E-4
    N_coeff = np.zeros(len(sn_range))

    path = '/home/xzk/work/github/MW-LMC-SCF/code/data_KL_SN/'
    S, T = coefficients_smoothing.read_coeff_matrix(path+'mwlmc_hal_sn_test_coeff_sample_', i+1, 20, 20, 20, 0, i+1, snaps=1)
    SS, TT, ST = coefficients_smoothing.read_cov_elements(path+'mwlmc_hal_sn_test_covmat_sample_', i+1, 20, 20, 20, 0, i+1, snaps=1)
    S_0, T_0, N_smooth_0 = coefficients_smoothing.smooth_coeff_matrix(S, T, SS, TT, ST, m_p_sim, 20, 20, 20, 0)# sn
    xyz1 = np.array([pos[:i*(particles_in_batch),0],
                    pos[:i*(particles_in_batch),1],
                    pos[:i*(particles_in_batch),2]]).T

    xyz2 = np.array([pos[(i+1)*(particles_in_batch):,0],
                    pos[(i+1)*(particles_in_batch):,1],
                    pos[(i+1)*(particles_in_batch):,2]]).T

    xyz = np.concatenate((xyz1, xyz2))
    assert(len(xyz==990000))
    rho_reference = biff.density(np.ascontiguousarray(xyz.astype(float)), S_0, T_0, M=1, r_s=40.85)       
    print('Here')
    """
    for sn in range(len(sn_range)):

        S_smooth, T_smooth, N_coeff[sn] = coefficients_smoothing.smooth_coeff_matrix(S, T, SS, TT, ST, m_p_sim, 20, 20, 20, sn_range[sn])


        rho_estimate = biff.density(np.ascontiguousarray(xyz.astype(float)), S_smooth, T_smooth, M=1, r_s=40.85) 
        Hp[sn] = np.sum(m_p*np.log(rho_factor*np.abs(rho_estimate))) 
        Kpq[sn] = np.sum(m_p*np.log(rho_factor*np.abs(rho_reference))) - Hp[sn]
        Kpq_rho[sn] = np.sum(rho_factor*rho_reference*np.log(np.abs(rho_reference/rho_estimate))) 

    """
    return Hp, Kpq, Kpq_rho, N_coeff
    
if __name__ == "__main__":

    #test = kpq_sn_estimator(path, 1)
    num_cores=2
    output = Parallel(n_jobs=num_cores)( delayed(kpq_sn_estimator)(i) for i in range(2) )
    print(output)
    
