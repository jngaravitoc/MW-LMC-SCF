import numpy as np
import matplotlib.pyplot as plt
import biff
import coefficients_smoothing



y_grid = np.arange(-300, 300, 5.5)
z_grid = np.arange(-300, 300, 5.5)
y_grid, z_grid = np.meshgrid(y_grid, z_grid)
bins=110

xyz = np.ascontiguousarray(np.array([np.zeros(len(y_grid.flatten())),
                                     y_grid.flatten(), z_grid.flatten()]).T)

def acceleration_plot(S, T, SS, TT, ST, figname, cbar_name):
    sn = [0, 1, 2, 3, 4, 5]
    fig = plt.figure(figsize=(10, 14))
    for i in range(len(sn)):
            S_smooth, T_smooth, N_smooth = coefficients_smoothing.smooth_coeff_matrix(S, T, SS, TT, ST, mass, 20, 20, 20, sn[i])
            a_biff = biff.gradient(np.ascontiguousarray(xyz), S_smooth, T_smooth, M=1, r_s=40.85, G=1)
            a_biff_all = np.sqrt(a_biff[:,0]**2 + a_biff[:,1]**2 + a_biff[:,2]**2)
            
            a_biff_0 = biff.gradient(np.ascontiguousarray(xyz), np.array([[[S_smooth[0,0,0]], [0], [0]]]).T,
                                     np.array([[[T_smooth[0,0,0]],[0],[0]]]).T, M=1, r_s=40.85, G=1)
            
            a_biff_all_0 = np.sqrt(a_biff_0[:,0]**2 + a_biff_0[:,1]**2 + a_biff_0[:,2]**2)
            
            plt.subplot(3, 2, i+1)
            #levels = np.arange(-0.03, 0.12, 0.005)
            im = plt.contourf(y_grid, z_grid, ((a_biff_all/a_biff_all_0)-1).reshape(bins, bins), 40,
                             origin='lower', cmap='viridis')

            plt.text(-180, 160, 'Ncoeff={}'.format(N_smooth), color='w')
            plt.text(-180, 130, 'S/N={}'.format(sn[i]), color='w')

            plt.xlim(-200, 200)
            plt.ylim(-200, 200)
            if ((i==4) | (i==5)):
                plt.xlabel('y[kpc]')
            if (i+1)%2==1:
                plt.ylabel('z[kpc]')
    cb_ax = fig.add_axes([0.93, 0.1, 0.02, 0.8])
    cbar = plt.colorbar(im, cax=cb_ax)
    fig.suptitle('Relative acceleration of MW + LMC unbound particles {}'.format(cbar_name), y=0.93)
    cbar.set_label('$\Delta |a|$')
    plt.savefig(figname+'.pdf', bbox_inches='tight')
    plt.savefig(figname+'.png', bbox_inches='tight')
    return 0
    
def relative_pot(S, T, SS, TT, ST, figname, cbar_name):
    sn = [0, 1, 2, 3, 4, 5]
    fig = plt.figure(figsize=(10, 14))
    for i in range(len(sn)):
            S_smooth, T_smooth, N_smooth = coefficients_smoothing.smooth_coeff_matrix(S, T, SS, TT, ST, mass, 20, 20, 20, sn[i])
            pot_biff = biff.potential(np.ascontiguousarray(xyz), S_smooth, T_smooth, M=1, r_s=40.85, G=1)
            pot_biff_0 = biff.potential(np.ascontiguousarray(xyz), np.array([[[S_smooth[0,0,0]], [0], [0]]]).T,
                                       np.array([[[T_smooth[0,0,0]],[0],[0]]]).T, M=1, r_s=40.85, G=1)


            plt.subplot(3, 2, i+1)
            #levels = np.arange(-4, -1, 0.1)
            im = plt.contourf(y_grid, z_grid, ((pot_biff/pot_biff_0)-1).reshape(bins, bins), 40,
                              origin='lower', cmap='inferno')

            plt.text(-180, 160, 'Ncoeff={}'.format(N_smooth), color='w')
            plt.text(-180, 130, 'S/N={}'.format(sn[i]), color='w')

            plt.xlim(-200, 200)
            plt.ylim(-200, 200)
            if ((i==4) | (i==5)):
                plt.xlabel('y[kpc]')
            if (i+1)%2==1:
                plt.ylabel('z[kpc]')

    cb_ax = fig.add_axes([0.93, 0.1, 0.02, 0.8])
    cbar = plt.colorbar(im, cax=cb_ax)
    fig.suptitle('Relative potential of MW + LMC unbound particles {}'.format(cbar_name), y=0.93)
    cbar.set_label('$\Delta \Phi$')
    plt.savefig(figname+'.pdf', bbox_inches='tight')    
    plt.savefig(figname+'.png', bbox_inches='tight')
    return 0
    
def rho_plots(S, T, SS, TT, ST, figname, cbar_name):
    sn = [0, 1, 2, 3, 4, 5]
    fig = plt.figure(figsize=(10, 14))
    for i in range(len(sn)):
            S_smooth, T_smooth, N_smooth = coefficients_smoothing.smooth_coeff_matrix(S, T, SS, TT, ST, mass, 20, 20, 20, sn[i])
            rho_biff = biff.density(np.ascontiguousarray(xyz), S_smooth, T_smooth, M=1, r_s=40.85)
            plt.subplot(3, 2, i+1)
            #levels = np.arange(-14, -2.4, 0.5)
            im = plt.contourf(y_grid, z_grid, np.log10(np.abs(rho_biff).reshape(bins, bins)), 40,
                              origin='lower', cmap='Spectral_r', )

            plt.text(-180, 160, 'Ncoeff={}'.format(N_smooth), color='k')
            plt.text(-180, 130, 'S/N={}'.format(sn[i]), color='k')

            plt.xlim(-200, 200)
            plt.ylim(-200, 200)
            if ((i==4) | (i==5)):
                plt.xlabel('y[kpc]')
            if (i+1)%2==1:
                plt.ylabel('z[kpc]')
    cb_ax = fig.add_axes([0.93, 0.1, 0.02, 0.8])
    cbar = plt.colorbar(im, cax=cb_ax)
    fig.suptitle('Density of MW + LMC unbound particles {}'.format(cbar_name), y=0.93)
    cbar.set_label(r'$\rm{Log_{10}} \rho$')
    plt.savefig(figname + '.pdf', bbox_inches='tight')
    plt.savefig(figname + '.png', bbox_inches='tight')
    return 0
    
def relative_rho(S, T, SS, TT, ST, figname, title_name):
    sn = [0, 1, 2, 3, 4, 5]
    fig = figure(figsize=(10, 14))
    for i in range(len(sn)):
            S_smooth, T_smooth, N_smooth = coefficients_smoothing.smooth_coeff_matrix(S, T, SS, TT, ST, mass, 20, 20, 20, sn[i])
            rho_biff = biff.density(np.ascontiguousarray(xyz), S_smooth, T_smooth, M=1, r_s=40.85)
            print(S_smooth[0,0,0])
            rho_biff_0 = biff.density(np.ascontiguousarray(xyz), np.array([[[S_smooth[0,0,0]], [0], [0]]]).T,
                                      np.array([[[T_smooth[0,0,0]],[0],[0]]]).T, M=1, r_s=40.85)
            subplot(3, 2, i+1)
            #levels = np.arange(-14, -2.4, 0.5)
            im = contourf(y_grid, z_grid, ((rho_biff/rho_biff_0) - 1).reshape(bins, bins), 30,
                          origin='lower', cmap='Spectral_r' )

            text(-180, 160, 'Ncoeff={}'.format(N_smooth), color='k')
            text(-180, 130, 'S/N={}'.format(sn[i]), color='k')

            xlim(-200, 200)
            ylim(-200, 200)
            if ((i==4) | (i==5)):
                xlabel('y[kpc]')
            if (i+1)%2==1:
                ylabel('z[kpc]')
                
    cb_ax = fig.add_axes([0.93, 0.1, 0.02, 0.8])
    cbar = colorbar(im, cax=cb_ax)
    fig.suptitle('Relative density of MW + LMC unbound particles {}'.format(title_name), y=0.93)
    cbar.set_label(r'$\Delta \rho$')
    savefig(figname + '.pdf', bbox_inches='tight')
    savefig(figname + '.png', bbox_inches='tight')    

if __name__ == "__main__":


    S_2T_V, T_2T_V = coefficients_smoothing.read_coeff_matrix('../data/MW/MW_lmc_unbound/mwlmc_unbound_BFE_2T_V_1e6_300_coeff_sample_0', 40, 20, 20, 20, 0, 40)
    SS_2T_V, TT_2T_V, ST_2T_V = coefficients_smoothing.read_cov_elements('../data/MW/MW_lmc_unbound/mwlmc_unbound_BFE_2T_V_1e6_300_covmat_sample_0', 40, 20, 20, 20, 0, 40)


    S_T_V, T_T_V = coefficients_smoothing.read_coeff_matrix('../data/MW/MW_lmc_unbound/mwlmc_unbound_BFE_T_V_1e6_300_coeff_sample_0', 40, 20, 20, 20, 0, 40)
    SS_T_V, TT_T_V, ST_T_V = coefficients_smoothing.read_cov_elements('../data/MW/MW_lmc_unbound/mwlmc_unbound_BFE_T_V_1e6_300_covmat_sample_0', 40, 20, 20, 20, 0, 40)

    S_T_2V, T_T_2V = coefficients_smoothing.read_coeff_matrix('../data/MW/MW_lmc_unbound/mwlmc_unbound_BFE_T_2V_1e6_300_coeff_sample_0', 40, 20, 20, 20, 0, 40)
    SS_T_2V, TT_T_2V, ST_T_2V = coefficients_smoothing.read_cov_elements('../data/MW/MW_lmc_unbound/mwlmc_unbound_BFE_T_2V_1e6_300_covmat_sample_0', 40, 20, 20, 20, 0, 40)

    S_T_4V, T_T_4V = coefficients_smoothing.read_coeff_matrix('../data/MW/MW_lmc_unbound/mwlmc_unbound_BFE_T_4V_1e6_300_coeff_sample_0', 40, 20, 20, 20, 0, 40)
    SS_T_4V, TT_T_4V, ST_T_4V = coefficients_smoothing.read_cov_elements('../data/MW/MW_lmc_unbound/mwlmc_unbound_BFE_T_4V_1e6_300_covmat_sample_0', 40, 20, 20, 20, 0, 40)


    mass = 1.577212515257997438e-06
    
    #acceleration_plot(S_2T_V, T_2T_V, SS_2T_V, TT_2T_V, ST_2T_V, 'delta_a_mwlmc_unbound_2T_V_sn' , '2T$<$V')
    #acceleration_plot(S_T_V, T_T_V, SS_T_V, TT_T_V, ST_T_V, 'delta_a_mwlmc_unbound_T_V_sn' , 'T$<$V')
    #relative_pot(S_T_2V, T_T_2V, SS_T_2V, TT_T_2V, ST_T_2V, 'delta_pot_mwlmc_unbound_T_2V_sn' , 'T$<$2V')
    #relative_pot(S_T_4V, T_T_4V, SS_T_4V, TT_T_4V, ST_T_4V, 'delta_pot_mwlmc_unbound_T_4V_sn' , 'T$<$4V')

    relative_rho(S_2T_V, T_2T_V, SS_2T_V, TT_2T_V, ST_2T_V, 'delta_a_mwlmc_unbound_2T_V_sn' , '2T$<$V')
    relative_rho(S_T_V, T_T_V, SS_T_V, TT_T_V, ST_T_V, 'delta_a_mwlmc_unbound_T_V_sn' , 'T$<$V')
    relative_rho(S_T_2V, T_T_2V, SS_T_2V, TT_T_2V, ST_T_2V, 'delta_pot_mwlmc_unbound_T_2V_sn' , 'T$<$2V')
    relative_rho(S_T_4V, T_T_4V, SS_T_4V, TT_T_4V, ST_T_4V, 'delta_pot_mwlmc_unbound_T_4V_sn' , 'T$<$4V')

