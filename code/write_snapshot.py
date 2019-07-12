"""
Usage:
    write_snapshot.py filename(has to be an .hdf5 extention), nparticles
"""


import numpy as np
import pygadgetreader
import gadget
import sys

filename = sys.argv[1]
particle_data = sys.argv[2]


def write_ics(filename, nparticles, pos, vel, mass):
   ics = gadget.ICs(filename, [0,nparticles], verbose=True)
   ics.pos[:] = pos
   ics.vel[:] = vel
   ics.mass[:] = np.ones(nparticles)*mass
   ics.ParticleIDs[:]=np.arange(1, nparticles+1)
   ics.write()


data = np.loadtxt(particle_data)


pos = np.array([data[:,0], data[:,1], data[:,2]]).T
vel = np.array([data[:,3], data[:,4], data[:,5]]).T
mass = data[0,6]


write_ics(filename, len(pos), pos, vel, mass)
