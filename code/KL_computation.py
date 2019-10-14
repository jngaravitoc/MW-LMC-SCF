import numpy as np
from joblib import Parallel, delayed



def write_output(file_name, data):
    np.savetxt(file_name, data)
	
def DKL(rho_i, rho_j):
    Dkl_ij = np.sum(mass * np.log10(np.abs(rho_i/rho_j)))
    return Dkl_ij 
	
def KL_batches(i, sn):
    Dkl_i = np.zeros(nbatches-1)
    rho = np.loadtxt(densities_file+'_{:0>3d}_sn_{:0>3d}.txt'.format(i, sn))
    k = 0
    for j in range(nbatches):
        if i!=j:
            Dkl_i[k] = DKL(rho[:,i], rho[:,j])
            k+=1
    print(Dkl_i)
    return np.sum(Dkl_i)/(nbatches-1) 

def DKL_all(sn):
    Dkl = np.zeros(nbatches)
    for i in range(nbatches):
        Dkl[i] = KL_batches(i, sn)
    
    print(np.sum(Dkl)/nbatches, sn)
    return np.sum(Dkl)/nbatches
	


if __name__ == "__main__":
    densities_file = "./data_KL_SN/rho_estimates/MW_100M_b1_dm_part_1e6_300_KL_analysis_batch"
    file_name = "DKL_sn.txt"
    npart = 10000
    mass = 1 / npart
    nbatches = 4
    sn_range = np.arange(0, 8, 0.2)
    num_cores=4


    output = Parallel(n_jobs=num_cores)(delayed(DKL_all)(sn) for sn in range(len(sn_range)))
    print(output)
    write_output(file_name, output)

