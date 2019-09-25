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

	contours = np.linspace(np.nanmin(np.log10(np.abs(dens))), np.nanmax(np.log10(np.abs(dens))), nbins)
	index_dens1 = []
	r_shell_mean = []
	N_dots_r = []

	for i in range(1, len(contours)-1):
		delta_rho_low = (contours[i+1]-contours[i])/2.
		delta_rho_high = (contours[i]-contours[i-1])/2.
		index_dens =  np.where((np.log10(np.abs(dens))>=contours[i] - delta_rho_low)
		                       & (np.log10(np.abs(dens))<=contours[i] + delta_rho_high))
		
		r_shell = (x_grid[index_dens]**2 + y_grid[index_dens]**2 + z_grid[index_dens]**2)**0.5
		max_r =  np.max(r_shell)

		N_dots = len(index_dens[0])
		if (max_r < rmax) & (N_dots > N_min):
		    
		    index_dens1.append(index_dens)
		    r_shell_mean.append(np.median(r_shell))
		    N_dots_r.append(N_dots)

	return index_dens1, r_shell_mean, N_dots_r


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
	axis = (3*eigval/N_dots)**0.5
	return eigvec, axis, s, q


#def plot_halo_ellipsoid():


#def goodness_of_fit():

