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
import datetime
import coefficients_smoothing

## TODO : work with gala
## Load S and T

def load_scf_coefficients(coeff_files, cov_files, nmax, lmax, 
        mmax, min_sample, max_sample, pmass, sn):
    """
    Load coefficients.
    TODO : REMOVE THIS FUNCTION WHEN SCRIPT IS FULLY WORKING 
    """
    
    ncoeff_sample = max_sample - min_sample 
    print(pmass)
    S, T = coefficients_smoothing.read_coeff_matrix(
            coeff_files, ncoeff_sample, nmax, lmax, mmax,
            min_sample, max_sample, snaps=90)

    SS, TT, ST = coefficients_smoothing.read_cov_elements(
            cov_files, ncoeff_sample, nmax, lmax, mmax, 
            min_sample, max_sample, snaps=90)
    
    S_smooth, T_smooth, N_smooth = coefficients_smoothing.smooth_coeff_matrix(
            S, T, SS, TT, ST, pmass, nmax, lmax, mmax, sn, sn_out=0)

    return S_smooth, T_smooth


def load_scf_coefficients_cov(coeff_files, nmax, lmax, 
        mmax, min_sample, max_sample, pmass, sn):
    """
    Load smoothed coefficients

    Parameters:
    =========
    coeff_files : str
        path to coefficients files
    nmax : int
        nmax order of the expansion
    lmax : int
        lmax order of the expansion
    mmax : int   
        mmax order of the expansion
    min_sample : int
        minimum value of variance sample 
    max_sample : int
        maximum value of variance sample 
    pmass : float
        particle mass
    sn : float
        signal-to-noise

    Return:
    =======

    Ssmooth : smoothed coefficients
    Tsmooth : smoothed coefficients

    """
    
    ncoeff_sample = max_sample - min_sample 
    print(pmass)
    S, T, SS, TT, ST = coefficients_smoothing.read_coeffcov_matrix(
            coeff_files, ncoeff_sample, nmax, lmax, mmax,
            min_sample, max_sample, snaps=0)
    
    S_smooth, T_smooth, N_smooth = coefficients_smoothing.smooth_coeff_matrix(
            S, T, SS, TT, ST, pmass, nmax, lmax, mmax, sn)

    return S_smooth, T_smooth


## grid to compute densities
def grid(box_size, ngrid_points, dim):
    r_grid = np.linspace(-box_size/2, box_size/2, ngrid_points)
    if dim==2:
        x_grid, y_grid = np.meshgrid(r_grid, r_grid)
        nbins = len(y_grid)
        return x_grid, y_grid, nbins
    elif dim==3:
        x_grid, y_grid, z_grid = np.meshgrid(r_grid, r_grid, r_grid)
        nbins = len(y_grid)
        return x_grid, y_grid, z_grid

def grid_spherical(N, rin, rout):
    ncomp = int(N**(1/3.))
    print(ncomp)
    phi = 2*np.linspace(0, 1, ncomp)*np.pi
    theta = np.arccos(2*np.linspace(0, 1, ncomp)-1)
    r = rin + np.linspace(0, 1, ncomp)**1/3.  * (rout-rin)
    R, T, P = np.meshgrid(r, theta, phi)
    x = R*np.sin(T)*np.sin(P)
    y = R*np.sin(T)*np.cos(P)
    z = R*np.cos(T)
    return np.ascontiguousarray(np.array([x.flatten(), y.flatten(), z.flatten()]).T)


## Compute densities
def scf_density_grid(
        x, y, z, Smw, Tmw, Slmc, Tlmc, nbins, rs_mw, rs_lmc, 
        lmc_com, quantity, G=1):
    
    q_all = np.zeros((nbins, nbins, nbins))
    for i in range(nbins):
        for j in range(nbins):
            for k in range(nbins):
                # These line is flipping x with y!
                xyz_lmc = np.array(
                        [[x[0,i,0]-lmc_com[0]], 
                        [y[j,0,0]-lmc_com[1]], 
                        [z[0,0,k]-lmc_com[2]]]).T
                xyz = np.array(
                        [[x[0,i,0]], 
                        [y[j,0,0]], 
                        [z[0,0,k]]]).T

                if quantity == "density":
                    q_all[i][j][k] = \
                        biff.density(xyz_lmc, Slmc, Tlmc, M=1, r_s=rs_lmc) \
                        + biff.density(xyz, Smw, Tmw, M=1, r_s=rs_mw)

                if quantity == "potential":
                    q_all[i][j][k] = \
                        biff.potential(np.ascontiguousarray(xyz_lmc), Slmc, Tlmc, M=1, r_s=rs_lmc,
                                G=G) \
                        + biff.potential(np.ascontiguousarray(xyz), Smw, Tmw, M=1, r_s=rs_mw, G=G)

                if quantity == "acceleration":
                    almc = biff.gradient(xyz_lmc, Slmc, Tlmc, M=1, r_s=rs_lmc, G=G) 
                    amw = biff.gradient(xyz, Smw, Tmw, M=1, r_s=rs_mw, G=G)
                    q_all[i][j][k] = np.sqrt(np.sum(almc**2)) \
                        + np.sqrt(np.sum(amw**2))

    return q_all.flatten()

def scf_density_grid_fast(
        x, y, z, Smw, Tmw, Slmc, Tlmc, nbins, rs_mw, rs_lmc, 
        lmc_com, quantity, G=1):
    
    q_all = np.zeros((nbins, nbins, nbins))
    xyz_lmc = np.array([x-lmc_com[0], 
                        y-lmc_com[1], 
                        z-lmc_com[2]]).T
    xyz = np.array([x,y,z]).T
    print(len(xyz_lmc[:,0]))
    print(xyz_lmc[:,0])
    if quantity == "density":
        q_all = biff.density(np.ascontiguousarray(xyz_lmc), Slmc, Tlmc, M=1, r_s=rs_lmc) \
                + biff.density(np.ascontiguousarray(xyz), Smw, Tmw, M=1, r_s=rs_mw)

    if quantity == "potential":
        q_all = biff.potential(np.ascontiguousarray(xyz_lmc), Slmc, Tlmc, M=1, r_s=rs_lmc, G=G) \
                + biff.potential(np.ascontiguousarray(xyz), Smw, Tmw, M=1, r_s=rs_mw, G=G)

    if quantity == "acceleration":
        almc = biff.gradient(xyz_lmc, Slmc, Tlmc, M=1, r_s=rs_lmc, G=G) 
        amw = biff.gradient(xyz, Smw, Tmw, M=1, r_s=rs_mw, G=G)
        q_all = np.sqrt(np.sum(almc**2, axis=1)) + np.sqrt(np.sum(amw**2, axis=1))

    return q_all


## Compute densities general 
## TODO: add off-set
def scf_density(x_grid, y_grid, z_grid, S, T, r_s_mw):
    xyz = np.ascontiguousarray(np.double(
        np.array([x_grid.flatten(), y_grid.flatten(), z_grid.flatten()]).T))

    dens_ratio_all = biff.density(xyz, S, T, M=1, r_s=r_s_mw)

    return dens_ratio_all

def scf_potential(x_grid, y_grid, z_grid, S, T, r_s_mw):
    xyz = np.ascontiguousarray(np.double(
        np.array([x_grid.flatten(), y_grid.flatten(), z_grid.flatten()]).T))

    pot_ratio_all = biff.potential(xyz, S, T, M=1, r_s=r_s_mw, G=1)

    return dens_ratio_all
## Write results.
def write_density(file_name, rho, x, y, z):    
    np.savetxt(file_name, np.array([rho, x, y, z]).T)
    return 0


if __name__ == "__main__":
    
    #scf_coef_files = '../data/interim/BFE/MWLMC3/bfe_MWLMC_unbound_81_MWLMC3_100M_new_b1_coeff_sample_'
    #scf_cov_files =  '../data/interim/BFE/MWLMC3/bfe_MWLMC_unbound_81_MWLMC3_100M_new_b1_covmat_sample_'
    
    #scf_coef_files_lmc = '../data/interim/BFE/MWLMC3/bfe_LMC3_bound_b0_coeff_sample_'
    #scf_cov_files_lmc = '../data/interim/BFE/MWLMC3/bfe_LMC3_bound_b0_covmat_sample_'
    
    #nmax = 20
    #lmax = 20
    #mmax = 20
    
    init = 0
    final = 10
    box_size = 400
    rmin = 20
    rmax = 30
    nbins = 512000
    rs_mw = 40.85
    rs_lmc = 10
    #pmass = 1.1996922E-6
    pmass = 1.577212515257997438e-06
    sn = 4
    sn_lmc = 10
    density_fname = "test_mwlmc_bfe_in_bs400_500_bins.txt"
    grid_dim = 3
    quantity = "density"
    G = 1

    xlmc_com = 7.57664148 - 9.19391488
    ylmc_com = 0.72792051 - 42.1269949
    zlmc_com = -31.24426422 - (-3.31300547)
    lmc_com = [xlmc_com, ylmc_com, zlmc_com]

     
    Smw, Tmw, N = coefficients_smoothing.get_coefficients('LMC5', 'isotropic', 'MW', snap=110, sn=4, mass=1.844E-6)
    #Slmc, Tlmc, N = coefficients_smoothing.get_coefficients('LMC5', 'isotropic', 'LMC', snap=110, sn=4, mass=1.2E-6)
    #Smwlmc, Tmwlmc, N = coefficients_smoothing.get_coefficients('LMC5', 'isotropic', 'MWLMC', snap=110, sn=4, mass=1.844E-6)
    Slmc = np.zeros((1,1,1))
    Tlmc = np.zeros((1,1,1))
    #Smono = np.zeros((1,1,1))
    #Tmono = np.zeros((1,1,1))
    #Smono[0,0,0] = Smw[0,0,0]
    #Tmono[0,0,0] = Tmw[0,0,0]
    #x, y, z = grid(500, 200, 3)
    #dens_all_fast = scf_density_grid_fast(x.flatten(), y.flatten(), z.flatten(), Smwlmc, Tmwlmc, Slmc, Tlmc, 200, rs_mw, rs_lmc, lmc_com, quantity, G)
    #pot_all_fast = scf_density_grid_fast(x.flatten(), y.flatten(), z.flatten(), Smwlmc, Tmwlmc, Slmc, Tlmc, 200, rs_mw, rs_lmc, lmc_com, 'potential', G)
    #density_fname = "rho_mw5_bfe_{}_bs_{}.txt".format(nbins, 500)
    #write_density(density_fname, dens_all_fast, x.flatten(), y.flatten(), z.flatten())
    #density_fname = "pot_mw5_bfe_{}_bs_.txt".format(nbins, 500)
    #write_density(density_fname, pot_all_fast, x.flatten(), y.flatten(), z.flatten())
    
    for i in range(5, 301, 5):
        print("radius", i)
        pos = grid_spherical(nbins, i-2.5, i+2.5)
        x = pos[:,0]
        y = pos[:,1]
        z = pos[:,2]
        #dens = scf_density_mw(x, y, z, Ss, Ts, rs_mw)
    
        #print(print(datetime.datetime.now().time()))
        #dens_all = scf_density_grid(x, y, z, Ssmw, Tsmw, Sslmc, Tslmc, int(nbins**(1/3.))+1, rs_mw, rs_lmc, lmc_com, quantity, G)
        print(datetime.datetime.now().time())
        dens_all_fast = scf_density_grid_fast(x, y, z, Smw, Tmw, Slmc, Tlmc, int(nbins**(1/3.))+1, rs_mw, rs_lmc, lmc_com, quantity, G)
        pot_all_fast = scf_density_grid_fast(x, y, z, Smw, Tmw, Slmc, Tlmc,  int(nbins**(1/3.))+1, rs_mw, rs_lmc, lmc_com, 'potential', G)
        print(datetime.datetime.now().time())
        density_fname = "rho_mw5_bfe_{}_r_{}.txt".format(nbins, i)
        write_density(density_fname, dens_all_fast, x, y, z)
        density_fname = "pot_mw5_bfe_{}_r_{}.txt".format(nbins, i)
        write_density(density_fname, pot_all_fast, x, y, z)
    
