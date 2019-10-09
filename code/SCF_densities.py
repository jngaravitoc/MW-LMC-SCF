"""
Script to compute BFE densities, potentials and accelerations
in a grid using the coefficients of the BFE.

Written by Nicolas Garavito-Camargo
email: jngaravitoc@email.arizona.edu
08/10/19


Things to do:
------------

1. Compute the accelerations
2. Include the potential/density/acceleration of the disk.
3. Generalize the grid definitions.


"""

import numpy as np
import biff
import coefficients_smoothing


## Load S and T
def load_scf_coefficients(coeff_files, cov_files, nmax, lmax, mmax, min_sample, max_sample, pmass, sn):
    """
    Load coefficients.
    """
    
    ncoeff_sample = max_sample - min_sample 
    print(pmass)
    S, T = coefficients_smoothing.read_coeff_matrix(coeff_files, ncoeff_sample,\
                                                            nmax, lmax, mmax,\
                                                            min_sample, max_sample, snaps=0)
    SS, TT, ST = coefficients_smoothing.read_cov_elements(cov_files, ncoeff_sample,\
                                                            nmax, lmax, mmax, min_sample, max_sample, snaps=0)
    
    S_smooth, T_smooth, N_smooth = coefficients_smoothing.smooth_coeff_matrix(S, T, SS, TT, ST,\
                                                                              pmass, nmax, lmax, mmax, sn)
    return S_smooth, T_smooth

## grid to compute densities
def grid_density(box_size, ngrid_points):
    r_grid = np.linspace(-box_size/2, box_size/2, ngrid_points)
    x_grid, y_grid, z_grid = np.meshgrid(r_grid, r_grid, r_grid)
    nbins = len(r_grid)
    return x_grid, y_grid, z_grid



## Compute densities
def scf_density_mwlmc(x, y, z, Smw, Tmw, Slmc, Tlmc, nbins, rs_mw, rs_lmc, lmc_com):
    dens_all = np.zeros((nbins, nbins, nbins))
    for i in range(nbins):
        for j in range(nbins):
            for k in range(nbins):
                # These line is flipping x with y!
                dens_all[i][j][k] = biff.density(np.array([[x[0,i,0]-lmc_com[0]],\
                                                           [y[j,0,0]-lmc_com[1]],\
                                                           [z[0,0,k]-lmc_com[2]]]).T,
                                                 Slmc, Tlmc, M=1,r_s=rs_lmc) + \
                                    biff.density(np.array([[x[0,i,0]], [y[j,0,0]], [z[0,0,k]]]).T,
                                                 Smw, Tmw, M=1, r_s=rs_mw)
    return dens_all.flatten()


## Compute densities

def scf_density_mw(x_grid, y_grid, z_grid, S, T, r_s_mw):
    xyz = np.ascontiguousarray(np.double(np.array([x_grid.flatten(), y_grid.flatten(), z_grid.flatten()]).T))
    dens_ratio_all = biff.density(xyz, S, T, M=1, r_s=r_s_mw)
    return dens_ratio_all

## Write results.
def write_density(file_name, rho):
    
    np.savetxt(file_name, np.array(rho))
    return 0


if __name__ == "__main__":
    
    scf_coef_files = '../data/MW/MW_lmc_unbound/mwlmc_unbound_BFE_2T_V_1e6_300_coeff_sample_0'
    scf_cov_files =  '../data/MW/MW_lmc_unbound/mwlmc_unbound_BFE_2T_V_1e6_300_covmat_sample_0'
    
    scf_coef_files_lmc ='../data/LMC/BFE_bound/LMC_1M_bound_2T_V_BFE_coeff_sample_0'
    scf_cov_files_lmc = '../data/LMC/BFE_bound/LMC_1M_bound_2T_V_BFE_covmat_sample_0'
    
    nmax = 20
    lmax = 20
    mmax = 20
    
    init = 0
    final = 10
    box_size = 200
    nbins = 301
    rs_mw = 40.85
    rs_lmc = 10
    pmass = 1.1996922E-6
    sn = 4
    sn_lmc = 10
    density_fname = "rho_mwlmc_bfe_in_100.txt"
    
    xlmc_com = 7.57664148 - 9.19391488
    ylmc_com = 0.72792051 - 42.1269949
    zlmc_com = -31.24426422 - (-3.31300547)
    lmc_com = [xlmc_com, ylmc_com, zlmc_com]

    
    Ssmw, Tsmw = load_scf_coefficients(scf_coef_files, scf_cov_files, nmax, lmax, mmax, init, final, pmass, sn)
    Sslmc, Tslmc = load_scf_coefficients(scf_coef_files_lmc, scf_cov_files_lmc, nmax, lmax, mmax, init, final, pmass, sn_lmc)

    x, y, z = grid_density(box_size, nbins)
    #dens = scf_density_mw(x, y, z, Ss, Ts, rs_mw)
    dens_all = scf_density_mwlmc(x, y, z, Ssmw, Tsmw, Sslmc, Tslmc, nbins, rs_mw, rs_lmc, lmc_com)
    write_density(density_fname, dens_all)
