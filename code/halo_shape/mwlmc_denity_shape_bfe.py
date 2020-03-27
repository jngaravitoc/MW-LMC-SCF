"""
Script to compute the halo shape based on the density
of a density contour defined using the BFE.

"""
import numpy as np
import matplotlib.pyplot as plt
import gala
import jellyfish
import sys
sys.path.append('../')
import coefficients_smoothing


def load_mwlmc_coefficients(cov_path, coeff_path, mass, nfiles, nmax, lmax):
	"""
	Load mwlmc coefficients
	"""
	S_mwlmc, T_mwlmc, N_mwlmc = coefficients_smoothing.smooth_coeff(coeff_mwlmc_path, covmat_mwlmc_path, 
											0, nfiles, nmax, lmax, lmax, nfiles, mass)
	return S_mwlmc, T_mwlmc


def compute_bfe_density_grid(S, T, npoints_grid):
	#
	#rho = 
	#grid = 
	return 0

def load_density_grid(dens_file, nbins):
	rho = np.loadtxt(dens_file)
	rho_matrix = np.reshape(rho, (nbins, nbins, nbins))
	return rho_matrix


def grid_points(nbins, dmin, dmax, dim='3d'):
	"""
	"""
	x_grid = np.linspace(dmin, dmax, nbins)
	y_grid = np.linspace(dmin, dmax, nbins)
		
	if dim == '3d':
		z_grid = np.linspace(dmin, dmax, nbins)
		x_grid, y_grid, z_grid = np.meshgrid(x_grid, y_grid, z_grid)
		return x_grid, y_grid, z_grid

	if dim == '2d':
		x_grid, y_grid = np.meshgrid(x_grid, y_grid)
		return x_grid, y_grid


def compute_density_contour(dens, nbins, x_grid, y_grid, z_grid, rmax=300, N_min=300):
    """
    Compute density contours
    Parameters:
    -----------
    dens : numpy.ndarray
    density array. 
    nbins : int
    Number of density bins for the contours.
    x_grid : numpy.ndarray
    y_grid : numpy.ndarray
    z_grid : numpy.ndarray
    Returns:
    --------
    """

    assert np.shape(x_grid) == np.shape(y_grid) == np.shape(z_grid) == np.shape(dens)
    contours = np.linspace(np.nanmin(np.log10(np.abs(dens))), np.nanmax(np.log10(np.abs(dens))), nbins+2)
    index_dens1 = []
    r_shell_mean = []
    N_dots_r = []
    for i in range(1, len(contours)-2):
        delta_rho_low = (contours[i+1]-contours[i])/2.
        delta_rho_high = (contours[i]-contours[i-1])/2.
        index_dens = np.where((np.log10(np.abs(dens))>=contours[i] - delta_rho_low) & (np.log10(np.abs(dens))<=contours[i] + delta_rho_high))[0]
        
        N_dots = len(index_dens)
        if N_dots > N_min: 
            r_shell = (x_grid[index_dens]**2 + y_grid[index_dens]**2 + z_grid[index_dens]**2)**0.5
            max_r = np.max(r_shell)
        if (max_r < rmax) & (N_dots > N_min):
            index_dens1.append(index_dens)
            r_shell_mean.append(np.median(r_shell))
            N_dots_r.append(N_dots)
    return index_dens1, r_shell_mean, N_dots_r


def twod_fits_plot(xgrid, ygrid, zgrid, xgrid_fit, ygrid_fit, zgrid_fit):
    x_data1, y_data1 = jellyfish.shapes.twod_surface(xgrid, ygrid)
    x_data2, y_data2 = jellyfish.shapes.twod_surface(xgrid, zgrid)
    y_data3, z_data3 = jellyfish.shapes.twod_surface(ygrid, zgrid)
                 
    x_fit1, y_fit1 = jellyfish.shapes.twod_surface(xgrid_fit, ygrid_fit)
    x_fit2, y_fit2 = jellyfish.shapes.twod_surface(xgrid_fit, zgrid_fit)
    y_fit3, z_fit3 = jellyfish.shapes.twod_surface(ygrid_fit, zgrid_fit)
                                    
    return [x_data1, y_data1, x_data2, y_data2, y_data3, z_data3], [x_fit1, y_fit1, x_fit2, y_fit2, y_fit3, z_fit3]


def compute_halo_shape(dens_contour, x_grid, y_grid, z_grid):
    """
    Computes halo shape
    Parameters:
    -----------
    dens_contour : numpy.ndarray
    	Array with the indices of where are the density contours in the halo
    x_grid : numpy.ndarray
    y_grid : numpy.ndarray
    z_grid : numpy.ndarray
    Returns:
    --------
    eigvec : numpy.ndarray (3,3)
    	Eigen vectors
    axis_lengths : numpy.3darray():
    	Length of principal axis.
    """
    shell_pos = np.array([x_grid[dens_contour], y_grid[dens_contour], z_grid[dens_contour]]).T
    eigvec, eigval, s, q = jellyfish.axis_ratios(shell_pos)
    N_dots = len(shell_pos[:,0])
    print("Number of points in grid:", N_dots)
    axis = (3*eigval/N_dots)**0.5
    return eigvec, axis, s, q

def halo_shape_r(rho_c, x, y, z):
    N_contours = len(rho_c)
    eigves = np.zeros((N_contours, 3, 3))
    axis = np.zeros((N_contours, 3))
    s = np.zeros(N_contours)
    q = np.zeros(N_contours)

    pos_ell = []
    for i in range(N_contours):
        eigves[i], axis[i], s[i], q[i] = compute_halo_shape(rho_c[i], x, y, z)
        pos_elli = jellyfish.ellipse_3dcartesian(axis[i], eigves[i])
        pos_ell.append(pos_elli)

    return pos_ell, eigves, axis, s, q


def axis_ratios(r_ell, s, q, figname):
    fig, ax = plt.subplots(1, 2, figsize=(14, 5))
    ax[0].plot(r_ell, s, '-o', label='$s$')
    ax[0].plot(r_ell, q, '-o', label='$q$')
    ax[0].axhline(0.8, label=r'$s\ Vera-Ciro+13$', c='C0', ls=':')
    ax[0].axhline(0.9, label=r'$q\ Vera-Ciro+13$', c='C1', ls=':')

    #ax[0].set_ylim(0.6, 1.01)
    ax[0].set_xlim(0, 270)
    ax[0].legend()
    ax[0].set_xlabel(r'$r_{ell}\rm{[kpc]}$')
    ax[0].set_ylabel(r'$s, q$')


    ax[1].plot(r_ell, (1-q**2)/(1-s**2), '-o')
    ax[1].axhline(0.67, c='k', ls='--', alpha=0.5)
    ax[1].axhline(0.33, c='k', ls='--', alpha=0.5)
    ax[1].text(10, 0.05, r'$\rm{Oblate}$')
    ax[1].text(10, 0.36, r'$\rm{Triaxial}$')
    ax[1].text(10, 0.7, r'$\rm{Prolate}$')
    ax[1].axhline(0.52, c='C0', ls=':', lw=2)
    ax[1].text(150, 0.53, 'Vera-Ciro+13')

    ax[1].set_ylim(0, 1)
    ax[1].set_xlim(0, 270)
    ax[1].set_xlabel(r'$r_{ell}\rm{[kpc]}$')
    ax[1].set_ylabel(r'$T\ \rm{Triaxiality}$')

    plt.savefig(figname, bbox_inches='tight')
    plt.close()

def contour_plot(data_2d_all, figname):
    fig, ax = plt.subplots(1, 3, figsize=(20,6))

    for i in range(0,len(data_2d_all)):
        
            
        xd = data_2d_all[i][1]
        yd = data_2d_all[i][0]
             
        
        if i==0:
            ax[0].plot(data_2d_all[i][1], data_2d_all[i][0], lw=1, label=r'$\rm{Shape\ tensor}$', ls='-', c='k')

        elif i>0:
            ax[0].plot(data_2d_all[i][1], data_2d_all[i][0], lw=1, ls='-', c='k')

    smin = 300
    ax[0].set_xlim(-smin, smin)
    ax[0].set_ylim(-smin, smin)

    ax[0].legend(fontsize=20)


    for i in range(0, len(data_2d_all)):
        #ax[1].plot(data_2d_all[i][4], data_2d_all[i][5], lw=2, c='C0')
        ax[1].plot(data_2d_all[i][4], data_2d_all[i][5], lw=1, ls='-', c='k')
        
    ax[1].set_xlim(-smin, smin)
    ax[1].set_ylim(-smin, smin)    

    for i in range(0, len(data_2d_all)):
        #ax[2].plot(data_2d_all[i][3], data_2d_all[i][2], lw=2, c='C0')
        ax[2].plot(data_2d_all[i][3], data_2d_all[i][2], lw=1, ls='-', c='k')
       
    ax[2].set_xlim(-smin, smin)
    ax[2].set_ylim(-smin, smin)

    ax[2].set_xlabel(r'$y \rm{[kpc]}$')
    ax[2].set_ylabel(r'$z \rm{[kpc]}$')


    ax[1].set_xlabel(r'$x \rm{[kpc]}$')
    ax[1].set_ylabel(r'$z \rm{[kpc]}$')


    ax[0].set_ylabel(r'$y \rm{[kpc]}$')
    ax[0].set_xlabel(r'$x \rm{[kpc]}$') 
    plt.savefig(figname, bbox_inches='tight')
    plt.close()

def hausdorff_distance(xyz1, xyz2):
    dxy=0
    for i in range(len(xyz1)):
        d = np.sqrt((xyz1[i,0] - xyz2[:,0])**2
                    + (xyz1[i,1] - xyz2[:,1])**2 
                    + (xyz1[i,2] - xyz2[:,2])**2)
        
        dxy += np.min(d)    
    return dxy/len(xyz1)

def dMH(xyz1, xyz2):
    d1 = hausdorff_distance(xyz1, xyz2)
    d2 = hausdorff_distance(xyz2, xyz1)
    return max(d1, d2)   

def load_data(filename):
    """

    """
    data = np.loadtxt(filename)
    q_all = data[:,0]
    x = data[:,1]
    y = data[:,2]
    z = data[:,3]
    return q_all, x, y, z

def compute_contours(r_ell, contours, x, y, z, figname):
    R = []
    s = []
    q = []
    data_2d_all = []
    fit_2d_all = []
    for i in range(len(r_ell)):
        R.append(r_ell[i])
        eigves_rho, axis_rho, sr, qr = compute_halo_shape(contours[i], x, y, z)
        s.append(sr)
        q.append(qr)
        print('Done computing density halo shape')
        #index_maxp = np.argmax(N_dotsp)
        pos_ell = jellyfish.ellipse_3dcartesian(axis_rho, eigves_rho)
        data_2d, fit_2d = twod_fits_plot(x[contours[i]], y[contours[i]], z[contours[i]], pos_ell[:,:,0].flatten(), pos_ell[:,:,1].flatten(), pos_ell[:,:,2].flatten())
        data_2d_all.append(data_2d)
        fit_2d_all.append(fit_2d) 
    
    axis_ratios(np.array(R), np.array(s), np.array(q), "axis_ratios"+figname)    
    contour_plot(data_2d_all, 'contours_'+figname) 
    contour_plot(fit_2d_all, 'fit_contours_'+figname) 
    return 0

def compute_mhd(r_ell, contours, x, y, z):
    d_MHD_r = []
    R = []
    for i in range(len(r_ell)):
        R.append(r_ell[i])
        eigves, axis, sr, qr = compute_halo_shape(contours[i], x, y, z)
        print(sr, qr)
        pos_ell = jellyfish.ellipse_3dcartesian(axis, eigves)
        xyz_c = np.array([x[contours[i]], y[contours[i]], z[contours[i]]]).T
        xyz_f = np.array([pos_ell[:,:,0].flatten(), pos_ell[:,:,1].flatten(), pos_ell[:,:,2].flatten()]).T
        print("Computing MHD density")
        if len(xyz_c) > 10000:
            rand_pos = np.random.randint(0, len(xyz_c), 10000)
            dmhd = dMH(xyz_c[rand_pos], xyz_f)
        else :
            dmhd = dMH(xyz_c, xyz_f)
        d_MHD_r.append(dmhd)

    print(R, d_MHD_r)
    
    fig, ax = plt.subplots(1, 2, figsize=(10,6))
    ax[0].plot(R, d_MHD_r)
    ax[1].plot(np.array(R), np.array(d_MHD_r)/np.array(R))
    plt.savefig('rho_r_mhd_r.pdf', bbox_inches='tight')
    plt.close()
    return np.array(R), np.array(d_MHD_r)
                                          
if __name__ == "__main__":
    filename = sys.argv[1]    
    figname = sys.argv[2]    
    n_contours = 50
    rho, xrho, yrho, zrho = load_data(filename)
    rho_contours, r_ell, N_dots = compute_density_contour(rho, n_contours, xrho, yrho, zrho)
    compute_contours(r_ell, rho_contours, xrho, yrho, zrho, figname)
 
    """
    dens_filename = sys.argv[1]    
    pot_filename = sys.argv[2]    
    figname = 'figure_test.png'
    n_contours = 50


    s_rho = []
    q_rho = []
    R_rho = []
    s_pot = []
    q_pot = []
    R_pot = []
    data_2d_all = []
    fit_2d_all = []
    data_2d_allp = []
    fit_2d_allp = []
    
    d_MHD_r = []
    d_MHD_p = []

    #for i in range(5, 300, 5):
    rho_data = np.loadtxt(dens_filename+".txt")
    pot_data = np.loadtxt(pot_filename+".txt")
    rho_all = rho_data[:,0]
    xrho = rho_data[:,1]
    yrho = rho_data[:,2]
    zrho = rho_data[:,3]
    pot_all = pot_data[:,0]
    xpot = pot_data[:,1]
    ypot = pot_data[:,2]
    zpot = pot_data[:,3]

    print("Done loading data")
    rho_contours, r_ell, N_dots = compute_density_contour(rho_all, n_contours, xrho, yrho, zrho)
    pot_contours, r_ellp, N_dotsp = compute_density_contour(pot_all, n_contours, xpot, ypot, zpot)
    print(r_ell)
    print(len(r_ell))
    print(N_dots, N_dotsp)
    
    for i in range(len(r_ell)):
        R_rho.append(r_ell[i])
        eigves_rho, axis_rho, sr, qr = compute_halo_shape(rho_contours[i], xrho, yrho, zrho)
        s_rho.append(sr)
        q_rho.append(qr)
        print('Done computing density halo shape')
        #index_maxp = np.argmax(N_dotsp)
        R_pot.append(r_ellp[i])
        eigves_pot, axis_pot, sr_pot, qr_pot =  compute_halo_shape(pot_contours[i], xpot, ypot, zpot)
        s_pot.append(sr_pot)
        q_pot.append(qr_pot)

        print('Done computing potential halo shape')
        pos_ell = jellyfish.ellipse_3dcartesian(axis_rho, eigves_rho)
        #data_2d, fit_2d = twod_fits_plot(x[rho_contours[index_max]], y[rho_contours[index_max]], z[rho_contours[index_max]], pos_ell[:,:,0].flatten(), pos_ell[:,:,1].flatten(), pos_ell[:,:,2].flatten())
        #data_2d_all.append(data_2d)
        #fit_2d_all.append(fit_2d)
    
        pos_ellp = jellyfish.ellipse_3dcartesian(axis_pot, eigves_pot)
        data_2dp, fit_2dp = twod_fits_plot(xp[pot_contours[index_maxp]], yp[pot_contours[index_maxp]], zp[pot_contours[index_maxp]], pos_ellp[:,:,0].flatten(), pos_ellp[:,:,1].flatten(), pos_ellp[:,:,2].flatten())
        #data_2d_allp.append(data_2dp)
        #fit_2d_allp.append(fit_2dp)
            
        xyz_c = np.array([xrho[rho_contours[i]], yrho[rho_contours[i]], zrho[rho_contours[i]]]).T
        xyz_f = np.array([pos_ell[:,:,0].flatten(), pos_ell[:,:,1].flatten(), pos_ell[:,:,2].flatten()]).T
        print(len(xyz_c))
        print(len(xyz_f))
        print("Computing MHD density")
        if len(xyz_c) > 10000:
            rand_pos = np.random.randint(0, len(xyz_c), 10000)
            dmhd = dMH(xyz_c[rand_pos], xyz_f)

        else :
            dmhd = dMH(xyz_c, xyz_f)
        d_MHD_r.append(dmhd)
        
                
        xyz_c = np.array([xpot[pot_contours[i]], ypot[pot_contours[i]], zpot[pot_contours[i]]]).T
        xyz_f = np.array([pos_ellp[:,:,0].flatten(), pos_ellp[:,:,1].flatten(), pos_ellp[:,:,2].flatten()]).T
        print("Computing MHD pot")
        if len(xyz_c) > 10000:
            rand_pos = np.random.randint(0, len(xyz_c), 10000)
            dmhd = dMH(xyz_c[rand_pos], xyz_f)
        else:
            dmhd = dMH(xyz_c, xyz_f)
        d_MHD_p.append(dmhd)

    
    plt.plot(R_rho, d_MHD_r)
    plt.savefig('rho_r_mhd.pdf', bbox_inches='tight')
    plt.close()
    plt.plot(np.array(R_rho), np.array(d_MHD_r)/np.array(R_rho))
    plt.savefig('rho_r_mhd_r.pdf', bbox_inches='tight')
    plt.close()
    
    plt.plot(R_rho, d_MHD_p)
    plt.savefig('pot_r_mhd.pdf', bbox_inches='tight')
    plt.close()
    plt.plot(np.array(R_rho), np.array(d_MHD_p)/np.array(R_rho))
    plt.savefig('pot_r_mhd_r.pdf', bbox_inches='tight')
    plt.close()
    
    #axis_ratios(np.array(R_rho), np.array(s_rho), np.array(q_rho), "rho"+figname)    
    #axis_ratios(np.array(R_pot), np.array(s_pot), np.array(q_pot), "pot"+figname)    
    contour_plot(data_2d_all, 'rho_contours_'+figname) 
    contour_plot(data_2d_allp, 'pot_contours_'+figname) 
    """
