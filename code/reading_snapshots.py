import numpy as np
from pygadgetreader import readsnap
from jellyfish import com


def host_sat_particles(xyz, vxyz, pids, pot, Nhost_particles, *p):
    """
    Function that return the host and the sat particles
    positions and velocities.

    Parameters:
    -----------
    xyz: snapshot coordinates with shape (n,3)
    vxys: snapshot velocities with shape (n,3)
    pids: particles ids
    Nhost_particles: Number of host particles in the snapshot
    Returns:
    --------
    xyz_mw, vxyz_mw, xyzlmc, vxyz_lmc: coordinates and velocities of
    the host and the sat.

    """
    sort_indexes = np.sort(pids)
    N_cut = sort_indexes[Nhost_particles]
    host_ids = np.where(pids<N_cut)[0]
    return xyz[host_ids], vxyz[host_ids], pids[host_ids], pot[host_ids], xyz[sat_ids], vxyz[sat_ids], pids[sat_ids], pot[sat_ids]

def host_particles(xyz, vxyz, pids, pot, mass, N_host_particles):
    """
    Function that return the host and the sat particles
    positions and velocities.

    Parameters:
    -----------
    xyz: snapshot coordinates with shape (n,3)
    vxys: snapshot velocities with shape (n,3)
    pids: particles ids
    Nhost_particles: Number of host particles in the snapshot
    
    Returns:
    --------
    xyz_mw, vxyz_mw, xyzlmc, vxyz_lmc: coordinates and velocities of
    the host and the sat.

    """
    sort_indexes = np.sort(pids)
    N_cut = sort_indexes[N_host_particles]
    host_ids = np.where(pids<N_cut)[0]
    return xyz[host_ids], vxyz[host_ids], pids[host_ids], pot[host_ids], mass[host_ids]


def sat_particles(xyz, vxyz, pids, pot, mass, Nhost_particles):
    """
    Function that return the host and the sat particles
    positions and velocities.

    Parameters:
    -----------
    xyz: snapshot coordinates with shape (n,3)
    vxys: snapshot velocities with shape (n,3)
    pids: particles ids
    Nhost_particles: Number of host particles in the snapshot
    Returns:
    --------
    xyz_mw, vxyz_mw, xyzlmc, vxyz_lmc: coordinates and velocities of
    the host and the sat.

    """
    sort_indexes = np.sort(pids)
    N_cut = sort_indexes[Nhost_particles]
    sat_ids = np.where(pids>=N_cut)[0]
    return xyz[sat_ids], vxyz[sat_ids], pids[sat_ids], pot[sat_ids], mass[host_ids]

def read_snap_coordinates(path, snap, N_halo_part, com_frame='MW', galaxy='MW'):
    """
    Returns the MW properties.
    
    Parameters:
    path : str
        Path to the simulations
    snap : name of the snapshot
    LMC : Boolean
        True or False if LMC is present on the snapshot.
    N_halo_part : int
        Number of particles in the MW halo.
    pot : boolean
        True or False if you want the potential back.
    com_frame : str
        Where the coordinates will be centered galactocentric (MW), on the
        satellite (sat), or in the LSR (LSR)
    galaxy : str
        galaxy coordinates to be returned (MW) or (sat)
    Returns:
    --------
    MWpos : 
    MWvel : 
    MWpot : 
    MWmass : 
    
    """
    # Load data
    all_pos = readsnap(path+snap, 'pos', 'dm')
    all_vel = readsnap(path+snap, 'vel', 'dm')
    all_ids = readsnap(path+snap, 'pid', 'dm')
    all_pot = readsnap(path+snap, 'pot', 'dm')


    
    
    print("Loading MW particles and LMC particles")

    if galaxy == 'MW':
        print('Loading MW particles')
        pos, vel, ids, pot, mass = host_particles(all_pos, all_vel, all_ids, all_pot, N_halo_part)

    elif galaxy == 'sat':
        print('Loading satellite particles')
        pos, vel, ids, pot, mass = sat_particles(all_pos, all_vel, all_ids, all_pot, N_halo_part)

    if com_frame == 'MW': 
        print('Computing coordinates in the MW COM frame')
        pos_disk = readsnap(path+snap, 'pos', 'disk')
        vel_disk = readsnap(path+snap, 'vel', 'disk')
        pot_disk = readsnap(path+snap, 'pot', 'disk')
        pos_cm, vel_cm = com.com_disk_potential(pos_disk, vel_disk, pot_disk)

    elif com_frame == 'sat':
        print('Computing coordinates in the satellite COM frame')
        if galaxy == 'MW':
            LMC_pos, LMC_vel, LMC_ids, LMC_pot, LMC_mass = sat_particles(all_pos, all_vel, all_ids, all_pot, all_mass, N_halo_part)
        pos_cm, vel_cm  = com.CM(LMC_pos, LMC_vel)

    elif com_frame=='LSR':
        print('Computing coordinates in the LSR frame')
        pos_disk = readsnap(path+snap, 'pos', 'disk')
        vel_disk = readsnap(path+snap, 'vel', 'disk')
        pot_disk = readsnap(path+snap, 'pot', 'disk')
        pos_cm, vel_cm = com.com_disk_potential(pos_disk, vel_disk, pot_disk)

        pos_LSR = np.array([-8.34, 0, 0])
        vel_LSR = np.array([11.1,  232.24,  7.25])
        
        print(pos_LSR)
        print(vel_LSR)
    assert len(MW_pos) == N_halo_part, 'something is wrong with the number of selected particles'

    return pos_cm, vel_cm, pot, mass
    




