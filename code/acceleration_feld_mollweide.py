import numpy as np
import matplotlib.pyplot as plt

# Librarires
import biff


# Astro
from astropy import units as u
import astropy.coordinates as coord
from astropy import constants

# Plotting
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from mpl_toolkits.axes_grid1.colorbar import colorbar
from mpl_toolkits.basemap import Basemap

# Routines 
# df 

import dynamical_friction as dyn_f
#import coefficients_smoothing
from mwlmc_bfe import smooth_coeff
import kinematics


class BFE_grid():
    #TODO: make COM positions an input parameter
    def __init__(self, distance):
        self.d = distance
        self.G_gadget = 43007.1
        self.rs = 40.85
        self.rs_sat = 10

    def spherical_grid_shift(self, xyz, shift_com):
        x_shift = xyz[:,0] - shift_com[0]
        y_shift = xyz[:,1] - shift_com[1]
        z_shift = xyz[:,2] - shift_com[2]
        return np.ascontiguousarray(np.array([x_shift, y_shift, z_shift]).T)

    def grid_cartessian(self):
        l4 = np.linspace(-179, 179, 200)
        b4 = np.linspace(-89, 89, 100)
        L4, B4 = np.meshgrid(l4, b4)
        Q4 = coord.Galactocentric(lon=L4*u.deg, lat=B4*u.deg, distance=self.d*u.kpc, representation='spherical')
        c1 = Q4.cartesian
        xyz = np.ascontiguousarray(np.array([c1.x.flatten(), c1.y.flatten(), c1.z.flatten()]).T)

        xyz_shift = self.spherical_grid_shift(xyz, [-1.6172733999999993, -41.39907439, -27.93125875])    
        return xyz, xyz_shift

    def a_bfe_all(self, S, T, S_mw, T_mw, S_mw_wake, T_mw_wake, S_lmc, T_lmc, S_mwlmc, T_mwlmc):
        xyz, xyz_shift = self.grid_cartessian()
        a_mon = biff.gradient(xyz, S, T, M=1, r_s=self.rs, G=self.G_gadget)
        a_mwwake = biff.gradient(xyz, S_mw, T_mw, G=self.G_gadget, M=1, r_s=self.rs)
        a_wake = biff.gradient(xyz, S_mw_wake, T_mw_wake, M=1, r_s=self.rs, G=self.G_gadget)
        a_mwlmc = biff.gradient(xyz, S_mwlmc, T_mwlmc, M=1, G=self.G_gadget, r_s=self.rs)
        # shift
        a_lmcshift = biff.gradient(xyz_shift, S_lmc, T_lmc, M=1, G=self.G_gadget, r_s=self.rs_sat)
        #a_lmcshift = biff.gradient(xyz, S_lmc, T_lmc, M=1, G=self.G_gadget, r_s=self.rs_sat)

        return a_mon, a_mwwake, a_wake, a_lmcshift, a_mwlmc

    def rho_bfe_all(self, S, T, S_mw, T_mw, S_mw_wake, T_mw_wake, S_lmc, T_lmc, S_mwlmc, T_mwlmc):
        xyz, xyz_shift = self.grid_cartessian()
        # Monopole
        rho_ref = biff.density(xyz, S_mw_ref, T_mw_ref, M=1, r_s=self.rs)        
        rho_wake = biff.density(xyz, S_mw, T_mw, M=1, r_s=self.rs)
        rho_mwwake = biff.density(xyz, S_mw_wake, T_mw_wake, M=1, r_s=self.rs)
        rho_lmcshift = biff.density(xyz_shift, S_lmc, T_lmc, M=1, r_s=self.rs_sat)
        #rho_lmcshift = biff.density(xyz, S_lmc, T_lmc, M=1, r_s=self.rs_sat)

        # MW + Debris
        rho_mwlmc = biff.density(xyz, S_mwlmc, T_mwlmc, M=1, r_s=self.rs)
        return rho_ref, rho_mwwake, rho_wake, rho_lmcshift, rho_mwlmc
        

def r_s(omegam, h, Mvir, c=False):
    rv = r_vir(omegam, h, Mvir)
    cv = c_vir(h, Mvir)
    if c:
        cv = c
    return rv/cv

# vdM 2012, equation A3
def f(x):
    a = np.log(1+x)
    b = x / (1+x)
    return a - b


def acceleration_LMC_analytic(M, c):
    self.rs = r_vir(0.3, 0.7, Mvir)/cvir
    self.x = r/self.rs
    x,y,z = position
    rr = np.sqrt(x**2. + y**2. + z**2.)

    return -G*self.Mvir*f(self.x)*i/(f(self.cvir)*self.rr**3.)


def a_hernquist(a, x, y, z, M, G):
    #G = constants.G
    #G = G.to(units.kiloparsec**3 / (units.Msun * units.s**2)) 
    G1 = G.to(u.kpc * u.m**2 / (u.Msun * u.s**2))
    print(G1)
    a = a * u.kpc
    x = x * u.kpc
    y = y * u.kpc
    z = z * u.kpc
    r = np.sqrt(x**2 + y**2 + z**2)
    M = M * u.Msun
    Ax = - 1.0 * x * G1 * M / (r * (r + a)**2)
    Ay = - 1.0 * y * G1 * M / (r * (r + a)**2)
    Az = - 1.0 * z * G1 * M / (r * (r + a)**2)
    print(min(Ax))
    Ax = Ax.to(u.km / u.s**2) 
    Ay = Ay.to(u.km / u.s**2)
    Az = Az.to(u.km / u.s**2)
    print(min(Ax))
    return [Ax, Ay, Az]

def acceleration_a_plot(ar, ap, at, atan, almc, amw, title='', v1=0, v2=0, fig_name=0, bar_label='', cbar_ticks=0):
    fig = plt.figure(figsize=(20, 5))
    plt.suptitle(title, y=0.85)
    ax1 = fig.add_subplot(132)
    m = Basemap(projection='moll',lon_0=0, lat_0=0, ax=ax1)
    ax1.set_title(r'$a_r$', fontsize=20)
    
    if ((v1!=0) & (v2!=0)):
        im = m.contourf(L4, B4, np.log10(np.abs(ar)).reshape(100, 200), 
                        200, extent=[-180, 180, -90, 90],
                        latlon=True, cmap='inferno_r', levels=np.linspace(v1, v2, 100), extend='both')
    else : 
        im = m.contourf(L4, B4, 
                        np.log10(np.abs(ar)).reshape(100, 200), 200, extent=[-180, 180, -90, 90],
                       latlon=True, cmap='inferno_r', extend='both')
    
    l_lmc, b_lmc = m(pos_lmc_gal[0]*180/np.pi, pos_lmc_gal[1]*180/np.pi)
    m.drawmeridians(np.arange(-180, 180, 180), linewidth=1.5, labels=[True, True, True], color='w')
    m.drawparallels(np.arange(-90, 90, 90), linewidth=1.5, color='w')
    m.scatter(l_lmc, b_lmc, marker='*', s=180, c='w', edgecolors='k', alpha=1) 
    
    ax2 = fig.add_subplot(131, projection='mollweide')
   
    ax2.set_title(r'$a_t$', fontsize=20)
    im2= ax2.streamplot(L4*np.pi/180, B4*np.pi/180, 
                       ap.reshape(100, 200), at.reshape(100, 200), 
                       color=(np.log10(atan).reshape(100, 200)),
                       linewidth=1, cmap=plt.cm.inferno_r)
    
    #xticks([])

    ax2.scatter(pos_lmc_gal[0], pos_lmc_gal[1], marker='*', s=180, c='k')
    ax2.set_xlabel(r'$\rm{Lon}[^{\circ}]$')
    ax2.set_ylabel(r'$\rm{Lat}[^{\circ}]$')

    
    ax1 = fig.add_subplot(133)
    m = Basemap(projection='moll',lon_0=0, lat_0=0, ax=ax1)
    ax1.set_title(r'$a_{lmc}/a_{mw}$', fontsize=20)

    im3 = m.contourf(L4, B4, np.log10(almc/amw).reshape(100,200), 
                    200, extent=[-180, 180, -90, 90],
                    latlon=True, cmap='inferno_r', 
                    levels=np.linspace(-1.3, 0.4, 20), extend='both')
    
    m.drawmeridians(np.arange(-180, 180, 180), linewidth=1.5, labels=[True, True, True], color='w')
    m.drawparallels(np.arange(-90, 90, 90), linewidth=1.5, color='w')
    
    cbar_ax1 = fig.add_axes([0.11, 0.1, 0.23, 0.05])
    cb1 = colorbar(im2.lines, cax=cbar_ax1, orientation='horizontal')
    
    cbar_ax2 = fig.add_axes([0.38, 0.1, 0.27, 0.05])
    cb2 = plt.colorbar(im, cax=cbar_ax2, orientation='horizontal')
    if ((v1!=0) & (v2!=0)):
        cb2.set_ticks(np.linspace(v1,v2,9))
    
    cbar_ax3 = fig.add_axes([0.67, 0.1, 0.27, 0.05])
    cb3 = plt.colorbar(im3, cax=cbar_ax3, orientation='horizontal')
    cb3.set_ticks([-1.2, -0.8, -0.4, 0, 0.4])
    

    if fig_name !=0:
        plt.savefig(fig_name, bbox_inches='tight')
        

#def load_coefficients(path, filename, sn):


if __name__ == "__main__":

    # BFE path
    covmat_lmc_path = '../data/interim/BFE/MWLMC3/bfe_LMC3_bound_b0_covmat_sample_'
    coeff_lmc_path = '../data/interim/BFE/MWLMC3/bfe_LMC3_bound_b0_coeff_sample_'

    covmat_mw_path = '../data/interim/BFE/MWLMC3/bfe_MW_81_MWLMC3_bound_b1_covmat_sample_'
    coeff_mw_path =  '../data/interim/BFE/MWLMC3/bfe_MW_81_MWLMC3_bound_b1_coeff_sample_'

    covmat_mwlmc_path = '../data/interim/BFE/MWLMC3/bfe_MWLMC_unbound_81_MWLMC3_100M_new_b1_covmat_sample_'
    coeff_mwlmc_path =  '../data/interim/BFE/MWLMC3/bfe_MWLMC_unbound_81_MWLMC3_100M_new_b1_coeff_sample_'

    # mass of particles in the BFE
    mass = 1.577212515257997438e-06

    # S/N  cuts
    S_mw, T_mw, N_mw = smooth_coeff(coeff_mw_path, covmat_mw_path, 0, 9, 20, 20, 20, 6, mass, snap=90, sn_out=0)
    S_lmc, T_lmc, N_lmc = smooth_coeff(coeff_lmc_path, covmat_lmc_path, 0, 9, 20, 20, 20, 6, mass, snap=90, sn_out=0)
    S_mwlmc, T_mwlmc, N_mwlmc = smooth_coeff(coeff_mwlmc_path, covmat_mwlmc_path, 0, 9, 20, 20, 20, 6, mass, snap=90, sn_out=0)

    # Units
    # from bfe are in in km/Gyr/s
    # this is the conversion factor to km/s^2
    gyr_to_s = 1*u.Gyr.to(u.s)
    acc_units_km_s2 = (1*u.kpc/u.Gyr**2).to(u.km/u.s**2)
    G_gadget = 43007.1 / 1E10 * u.km * u.kpc**2/u.Gyr/u.s/ u.Msun		
    print(gyr_to_s)

    # Grid
    l4 = np.linspace(-179, 179, 200)
    b4 = np.linspace(-89, 89, 100)
    L4, B4 = np.meshgrid(l4, b4)

    # Name conventions:
    #ref = monopole term
    # wake = BFE with all the terms but the monopole 
    # mwwake = MW DM halo
    # lmc = LMC DM halo truncated
    # a_mwdebris = MW DM halo + debris

    # Wake:
    S_mw_wake = np.copy(S_mw)
    T_mw_wake = np.copy(T_mw)
    S_mw_wake[0,0,0] = 0
    T_mw_wake[0,0,0] = 0

    # Monopole
    S_mw_ref = np.zeros(np.shape(S_mw))
    T_mw_ref = np.zeros(np.shape(T_mw))
    S_mw_ref[:,0,0] = S_mw[:,0,0]
    T_mw_ref[:,0,0] = T_mw[:,0,0]

    # Parameters analytic calculation
    r_s_sat = 10
    M_sat = 0.8E11
    dist = [20, 40, 60, 80, 100, 120, 140, 160]
    fig_name = 'MW_LMC_BFE_a_analytic.png'
    title = 'MW+LMC analytic [50 kpc]'
    analytic_lmc = 1
    wake = 0

    # LMC positions in galactocentric units

    pos_lmc_gal = np.array([-1.60984022, -0.59317509])

    for i in range(len(dist)):
        fig_name = 'MW_LMC_BFE_a_analytic_{:0>3d}.png'.format(dist[i])
        title = 'MW+LMC analytic [{:0>3d} kpc]'.format(dist[i])
        grid = BFE_grid(dist[i])
        xyz, xyz_shift = grid.grid_cartessian()
        Cuadrant4 = kinematics.Kinematics(xyz, xyz)
        a_ref, a_wake, a_mwwake, a_lmc, a_mwdebris = grid.a_bfe_all(S_mw_ref, T_mw_ref, 
                                                                    S_mw_wake, T_mw_wake, 
                                                                    S_mw, T_mw, S_lmc, T_lmc, S_mwlmc,
                                                                    T_mwlmc)
        #rho_ref50, rho_wake50, rho_mwwake50, rho_lmc50, rho_mwdebris50 = grid.rho_bfe_all(S_mw_ref, T_mw_ref,
        #			                                                                          S_mw_wake, T_mw_wake, 
        #		                                                                          S_mw, T_mw,
        #		                                                                          S_lmc, T_lmc, 
        #		                                                                          S_mwlmc, T_mwlmc)

        ar_ref, at_ref, ap_ref = Cuadrant4.vel_cartesian_to_galactic(xyz, a_ref)
        ar_lmc, at_lmc, ap_lmc = Cuadrant4.vel_cartesian_to_galactic(xyz, a_lmc)
        ar_mwdebris, at_mwdebris, ap_mwdebris = Cuadrant4.vel_cartesian_to_galactic(xyz, a_mwdebris)
        ar_mwwake, at_mwwake, ap_mwwake = Cuadrant4.vel_cartesian_to_galactic(xyz, a_mwwake)
        ar_wake, at_wake, ap_wake = Cuadrant4.vel_cartesian_to_galactic(xyz, a_wake)


        # Analytic calculation
        if analytic_lmc == 1:		
	        a_lmc = a_hernquist(r_s_sat, xyz_shift[:,0], xyz_shift[:,1], 
			          	        xyz_shift[:,2], M_sat, G_gadget)

	        ar_lmc, at_lmc, ap_lmc = Cuadrant4.vel_cartesian_to_galactic(xyz, -np.array(a_lmc).T)

	        a_tan_lmc = np.sqrt(at_lmc**2 + ap_lmc**2)
	        a_lmc_mag = np.sqrt(a_lmc[:,0]**2 + a_lmc[:,1]**2 + a_lmc[:,2]**2)

        elif analytic_lmc == 0:
	        a_tan_lmc = np.sqrt(at_lmc**2 + ap_lmc**2)/gyr_to_s
	        a_lmc_mag = np.sqrt(a_lmc[:,0]**2 + a_lmc[:,1]**2 + a_lmc[:,2]**2)/gyr_to_s
	        ar_lmc = ar_lmc / gyr_to_s
	        at_lmc = at_lmc / gyr_to_s
	        ap_lmc = ap_lmc / gyr_to_s

        a_tan_ref = np.sqrt(at_ref**2 + ar_ref**2)/gyr_to_s
        a_ref_mag = np.sqrt(a_ref[:,0]**2 + a_ref[:,1]**2 + a_ref[:,2]**2)/gyr_to_s

        if wake == 0:
            acceleration_a_plot(ar_lmc+ar_ref/gyr_to_s, -ap_lmc+(ap_ref/gyr_to_s), 
				                -at_lmc+(at_ref/gyr_to_s), 
				                a_tan_lmc+a_tan_ref,
	                            a_lmc_mag, a_ref_mag,
	                            v1=-14.4, v2=-13.4, 
				                title=title, 
				                fig_name=fig_name)
        elif wake == 1:
            acceleration_a_plot(ar_lmc+ar_ref/gyr_to_s, -ap_lmc+(ap_ref/gyr_to_s), 
				                -at_lmc+(at_ref/gyr_to_s), 
				                a_tan_lmc+a_tan_ref,
	                            a_lmc_mag, a_ref_mag,
	                            v1=-14, v2=-13.4, 
				                title=title, 
				                fig_name=fig_name)


