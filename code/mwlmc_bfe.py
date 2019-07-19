import numpy as np
import matplotlib.pyplot as plt
import biff
import coefficients_smoothing


def smooth_coeff(coeff_path, cov_path, ni, nf, nmax, lmax, mmax, sn, pmass):
    nfiles = nf-ni+1
    S, T = coefficients_smoothing.read_coeff_matrix(coeff_path,  nfiles, nmax, \
                                                   lmax, mmax, nmin=ni, nmax=nf)
    SS, TT, ST = coefficients_smoothing.read_cov_elements(cov_path,  nfiles, nmax,\
                                                         lmax, mmax, nmin=ni, nmax=nf)
    S_smooth, T_smooth, N_smooth = coefficients_smoothing.smooth_coeff_matrix(S, T, SS, TT, ST, pmass, nmax, lmax, mmax, sn)
    return S_smooth, T_smooth, N_smooth



def grid(box_size, nbins):
    y_grid = np.linspace(-box_size/2., box_size/2., nbins)
    z_grid = np.linspace(-box_size/2., box_size/2., nbins)
    nbins = len(y_grid)
    y_grid, z_grid = np.meshgrid(y_grid, z_grid)

    return y_grid, z_grid, nbins

def combine_bfe_rho(S1, T1, S2, T2, y_grid, z_grid, lmc_com, nbins):
    rho_mwlmc = np.zeros((nbins, nbins))
    for i in range(nbins):
        for j in range(nbins):
            rho_mwlmc[i][j] = biff.density(np.array([[0-lmc_com[0]],
                                           [y_grid[0][i]-lmc_com[1]],
                                           [z_grid[j,0]-lmc_com[2]]]).T,
                                           S2, T2, M=1,r_s=10) + \
                              biff.density(np.array([[0], [y_grid[0][i]], [z_grid[j,0]]]).T,
                                           S1, T1, M=1, r_s=40.85)
    
    return rho_mwlmc


def combine_bfe_pot(S1, T1, S2, T2, y_grid, z_grid, lmc_com, nbins):
    pot_mwlmc = np.zeros((nbins, nbins))
    for i in range(nbins):
        for j in range(nbins):
            pot_mwlmc[i][j] = biff.potential(np.array([[0-lmc_com[0]],
                                           [y_grid[0][i]-lmc_com[1]],
                                           [z_grid[j,0]-lmc_com[2]]]).T,
                                           S2, T2, M=1, r_s=10, G=1) + \
                              biff.potential(np.array([[0], [y_grid[0][i]], [z_grid[j,0]]]).T,
                                           S1, T1, M=1, r_s=40.85, G=1)
    
    return pot_mwlmc

def combine_bfe_a(S1, T1, S2, T2, y_grid, z_grid, lmc_com, nbins):
	a_mwlmc = np.zeros((nbins, nbins))
	for i in range(nbins):
		for j in range(nbins):
			a = biff.gradient(np.array([[0-lmc_com[0]],
                                       [y_grid[0][i]-lmc_com[1]],
                                       [z_grid[j,0]-lmc_com[2]]]).T,
							 S2, T2, M=1,r_s=10, G=1) + \
                biff.gradient(np.array([[0], [y_grid[0][i]], [z_grid[j,0]]]).T,
							 S1, T1, M=1, r_s=40.85, G=1)
			a_mwlmc[i][j] = (a[:,0]**2 + a[:,1]**2 + a[:,2]**2)**0.5	
    
	return a_mwlmc

def plot_bfe(figname, qbfe, y_grid, z_grid, cmap, cbar_label, Ncoeff, title):
	fig = plt.figure(figsize=(6, 5))
	plt.contourf(y_grid, z_grid, qbfe, 30, origin='lower', cmap=cmap)
	colorbar = plt.colorbar()
	colorbar.set_label(cbar_label)
	plt.xlabel(r'$\rm{y[kpc]}$')
	plt.ylabel(r'$\rm{z[kpc]}$')
	plt.title(title)
	#plt.text(-270, 260, 'Ncoeff:'+str(Ncoeff))
	plt.savefig(figname, bbox_inches='tight')
	return 0

#def disk_potential():

if __name__ == "__main__":
    covmat_lmc_path = '../data/LMC/BFE_bound/LMC_1M_iterative_bound_T_V_BFE_covmat_sample_0'
    coeff_lmc_path = '../data/LMC/BFE_bound/LMC_1M_iterative_bound_T_V_BFE_coeff_sample_0'
    covmat_mwlmc_path = '../data/MW/MW_lmc_unbound/mwlmc_unbound_BFE_T_V_iterative_1e6_300_covmat_sample_0'
    coeff_mwlmc_path = '../data/MW/MW_lmc_unbound/mwlmc_unbound_BFE_T_V_iterative_1e6_300_coeff_sample_0'

    sn_lmc  = 3
    ni_lmc = 0
    nf_lmc = 40

    ni_mw = 0
    nf_mw = 4
    sn_mw = 4.5
    mass = 1.7995383e-05
    S_lmc, T_lmc, N_lmc = smooth_coeff(coeff_lmc_path, covmat_lmc_path, ni_lmc,\
                                       nf_lmc, 20, 20, 20, sn_lmc, mass) 
    S_mw, T_mw, N_mw = smooth_coeff(coeff_mwlmc_path, covmat_mwlmc_path, ni_mw,\
                                       nf_mw, 20, 20, 20, sn_mw, mass) 
    y_grid, z_grid, nbins = grid(600, 118)

    xlmc_com = 7.57664148 - 9.19391488 
    ylmc_com = 0.72792051 - 42.1269949
    zlmc_com = -31.24426422 - (-3.31300547)
    lmc_com = [xlmc_com, ylmc_com, zlmc_com]

    #rho_mwlmc = combine_bfe_rho(S_mw, T_mw, S_lmc, T_lmc, y_grid, z_grid, lmc_com, nbins)
    pot_mwlmc = combine_bfe_pot(S_mw, T_mw, S_lmc, T_lmc, y_grid, z_grid, lmc_com, nbins)
    #a_mwlmc = combine_bfe_a(S_mw, T_mw, S_lmc, T_lmc, y_grid, z_grid, lmc_com, nbins)

    #plot_bfe('rho_bfe_mwlmc_300.pdf',np.log10(np.abs(rho_mwlmc.reshape(nbins, nbins))),
	#		 y_grid, z_grid, 'Spectral_r', r'$\rho$', N_lmc+N_mw, r'$\rm{Log_{10} \rho_{mw+lmc}}$')
    plot_bfe('pot_bfe_mwlmc_300_lmc_sn_5.pdf',np.log10(np.abs(pot_mwlmc.reshape(nbins, nbins))),
			 y_grid, z_grid, 'inferno', r'$\Phi$', N_lmc+N_mw, r'$\rm{Log_{10} \Phi_{mw+lmc}}$')
    #plot_bfe('a_bfe_mwlmc_300.pdf',np.log10(np.abs(a_mwlmc.reshape(nbins, nbins))),
 	#		 y_grid, z_grid, 'viridis', '$a$', N_lmc+N_mw, r'$\rm{Log_{10} a_{mw+lmc}}$')



    

    

    

    




