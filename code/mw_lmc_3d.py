import biff
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from mayavi import mlab


font = {'size':18, 'family':'serif'}
matplotlib.rc('font', **font)



def read_satellites_6d_coords(file_name):
    data = np.genfromtxt(file_name, dtype=None)
    vx = data['f7']
    vy = data['f9']
    vz = data['f11']
    x = data['f13']
    y = data['f15']
    z = data['f17']
    names = data['f0']
   
    v_mag = np.sqrt(vx**2 + vy**2 + vz**2)
    return x, y, z, vx/v_mag, vy/v_mag, vz/v_mag, names


def co_orbiting_sats(names):
    print('satellites co-orbiting around MW:')
    co_orb_satellites = [b'Hor1', b'Cra2', b'Car3', b'Ret2',
                         b'Dra1', b'Hyd1', b'Dra2', b'UMin1',
                         b'Fnx1', b'Leo1', b'Car1']

    index_name = np.zeros(len(co_orb_satellites))
    for i in range(len(co_orb_satellites)):
        index_name[i] = np.where(names==co_orb_satellites[i])[0]

    print(names[index_name.astype(int)])
    return index_name.astype(int)


def counter_orbiting_sats(names):
    print('satellites counter-orbiting around MW:')
    counter_orb_satellites = [b'Leo4', b'Aqu2', b'CanVen2', b'Scu1', b'Boo3']

    index_name = np.zeros(len(counter_orb_satellites))
    for i in range(len(counter_orb_satellites)):
        index_name[i] = np.where(names==counter_orb_satellites[i])[0]

    print(names[index_name.astype(int)])
    return index_name.astype(int)
   


if __name__ == "__main__":

    all_sats = read_satellites_6d_coords('../data/MW_satellites.dat')

    x_sats = all_sats[0]
    y_sats = all_sats[1]
    z_sats = all_sats[2]
    vx_sats = all_sats[3]
    vy_sats = all_sats[4]
    vz_sats = all_sats[5]
    names = all_sats[6]

    print(len(names))
    co_orb = co_orbiting_sats(names)
    counter_orb = counter_orbiting_sats(names)
    
    mlab.figure(bgcolor=(0.0, 0.0, 0.0))


    mlab.points3d(x_sats, y_sats, z_sats, color=(1,0,0),\
                  scale_factor=10, opacity=1)
    mlab.quiver3d(x_sats, y_sats, z_sats, vx_sats, vy_sats,\
                  vz_sats, color=(1,0,0), scale_factor=30, opacity=1)


    mlab.points3d(x_sats[co_orb], y_sats[co_orb], z_sats[co_orb],\
                 color=(1,1,0), scale_factor=10, opacity=1)
    mlab.quiver3d(x_sats[co_orb], y_sats[co_orb], z_sats[co_orb],\
                  vx_sats[co_orb], vy_sats[co_orb], vz_sats[co_orb],\
                  color=(1,1,0), scale_factor=30, opacity=1)

    mlab.points3d(x_sats[counter_orb], y_sats[counter_orb], z_sats[counter_orb],\
                 color=(1,0,1), scale_factor=10, opacity=1)
    mlab.quiver3d(x_sats[counter_orb], y_sats[counter_orb], z_sats[counter_orb],\
                  vx_sats[counter_orb], vy_sats[counter_orb], vz_sats[counter_orb],\
                  color=(1,0,1), scale_factor=30, opacity=1)


    mlab.points3d(0, 0, 0, 20,  mode='2dcircle', color=(1,1,1),\
                 scale_factor=1.0, opacity=1, resolution=100)

    mlab.savefig('mw_satellites_3d.png', size=(200, 200))

#mlab.close()
#Function that reads the N-body simulation orbit
"""
def reading_Nbody(snap_name):
    data = np.loadtxt(snap_name)
    #time = data[:,0]
    #Rgal = data[:,1]
    x_sat= data[:,6]
    y_sat = data[:,7]
    z_sat = data[:,8]
    x_gal = data[:,0]
    y_gal = data[:,1]
    z_gal = data[:,2]
    #Vgal = data[:,8]
    vx_sat = data[:,9]
    vy_sat = data[:,10]
    vz_sat = data[:,11]
    vx_gal = data[:,3]
    vy_gal = data[:,4]
    vz_gal = data[:,5]
    Rgal= np.sqrt((x_sat-x_gal)**2 + (y_sat-y_gal)**2 + (z_sat-z_gal)**2)
    Vgal= np.sqrt((vx_sat-vx_gal)**2 + (vy_sat-vy_gal)**2 + (vz_sat-vz_gal)**2)

    return Rgal, x_sat, y_sat, z_sat, x_gal, y_gal, z_gal, Vgal, vx_sat, vy_sat, vz_sat, vx_gal, vy_gal, vz_gal


r_s_sims = 40.85
xlmc_com = 1
ylmc_com = -41
zlmc_com = -28



coeff_c = np.loadtxt('../../SCF_tools/PCA/MWLMC5_coeff_20_20_100M_b1.txt')
S = coeff_c[:,0]
T = coeff_c[:,1]

S_matrix = np.zeros((21, 21, 21))
T_matrix = np.zeros((21, 21, 21))


counter = 0
for n in range(21):
    for l in range(21):
        for m in range(0, l+1):
            S_matrix[n][l][m] = S[counter]
            T_matrix[n][l][m] = T[counter]
            counter +=1


coeff_lmc = np.loadtxt('./LMC/coeff_rand_40lmc5_b1_1E6_1.txt')
S_lmc = coeff_lmc[:,0]
T_lmc = coeff_lmc[:,1]

S_lmc_1e6 = np.zeros((41, 21, 21))
T_lmc_1e6 = np.zeros((41, 21, 21))


counter = 0
for n in range(41):
    for l in range(21):
        for m in range(0, l+1):
            S_lmc_1e6[n][l][m] = S_lmc[counter]
            T_lmc_1e6[n][l][m] = T_lmc[counter]
            counter +=1


x_grid = np.arange(-300, 300, 5)
y_grid = np.arange(-300, 300, 5)
z_grid = np.arange(-300, 300, 5)
X_grid, Y_grid, Z_grid = np.meshgrid(x_grid, y_grid, z_grid)

"""

#a_ratio_all = np.zeros((110, 110))
#for i in range(110):
#    print(y_grid[i])
#    for j in range(110):
#        print(z_grid[j], zlmc_com)
#        xyz_lmc = (np.array([[0-xlmc_com], [y_grid[i]-ylmc_com],[z_grid[j]-zlmc_com]]))
#        print(xyz_lmc.astype(float))
#        print(biff.gradient(xyz_lmc.astype(float), S_lmc_1e6, T_lmc_1e6, M=11.41*1E6, r_s=10, G=1))

        #a_all_lmc = biff.gradient(np.array([[0-xlmc_com], [y_grid[0][i]-ylmc_com], [z_grid[j,0]-zlmc_com]]).T,
        #                      S_lmc_1e6, T_lmc_1e6, M=11.41*1E6,r_s=10, G=1) 
        #a_all_mw = biff.gradient(np.array([[0], [y_grid[0][i]], [z_grid[j,0]]]).T, S_matrix,
        #                      T_matrix, M=1, r_s=40.85, G=1)
                              
        #a_all[i][j] = np.sqrt(a_all_lmc[0][0]**2 + a_all_lmc[0][1]**2 +\
        #                      a_all_lmc[0][2]**2) +  np.sqrt(a_all_mw[0][0]**2 + \
        #                      a_all_mw[0][1]**2 + a_all_mw[0][2]**2)



#np.savetxt('a_mwlmc.txt', a_ratio_all)



#rho_mwlmc = np.loadtxt('rhomwlmc.txt')
#rho_mwlmc_base = np.loadtxt('rhomwlmc_base.txt')

#rho_matrix = np.reshape(rho_mwlmc, (120, 120, 120))
#rho_base = np.reshape(rho_mwlmc_base, (120, 120, 120))

### lmc orbit

"""
LMC5_b1 = '../../MW_anisotropy/data/orbits/LMC5_100Mb1_orbit.txt'
LMC5_b1_orbit = reading_Nbody(LMC5_b1)

x_sat = LMC5_b1_orbit[1]
y_sat = LMC5_b1_orbit[2]
z_sat = LMC5_b1_orbit[3]
x_gal = LMC5_b1_orbit[4]
y_gal = LMC5_b1_orbit[5]
z_gal = LMC5_b1_orbit[6]
"""


#mlab.figure(bgcolor=(0.0, 0.0, 0.0))

#mlab.plot3d(y_sat[:111]-y_gal[:111]+2.5, x_sat[:111]-x_gal[:111]+2.5, z_sat[:111]-z_gal[:111]+2.5, 
#            np.ones(len(x_sat[:111])), color=(1,0,0), line_width=200,
#            tube_radius=2, opacity=1)

#mlab.points3d(y_sat[111]-y_gal[111] +2.5, x_sat[111]-x_gal[111]+2.5, z_sat[111]-z_gal[111]+2.5,
#             100, color=(1,0,0), scale_factor=2, opacity=1)
#
#mlab.points3d(2.5, 2.5, 2.5,
#             100,  mode='2dcircle', color=(1,0,0), scale_factor=1.0, opacity=1,
#             resolution=40, line_width=180)

#tt = mlab.contour3d((rho_matrix/rho_base)-1,contours=12,
#                    opacity=.3, vmin=-0.45, vmax=0.45, extent=[-300, 300, -300, 300, -300, 300],
#                    transparent=True)#, colormap='viridis')



#for i in range(0, 1):
#    mlab.view(azimuth=i, elevation=90, distance=1200)
#    mlab.savefig('all_a.png',
#            size=(4000,4000))

#mlab.close()
