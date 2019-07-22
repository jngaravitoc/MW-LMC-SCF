import numpy as np
from joblib import Parallel, delayed


def write_output(file_name, data):
    np.savetxt(file_name, data)


def sum_rho_true(i, sn):
	
	rho = np.loadtxt(densities_file+'_{:0>3d}_sn_{:0>3d}.txt'.format(i, sn))
	
	rho_true = np.zeros(nbatches-1)
	k=0
	L1 = np.sum(rho[:,i]**2)
	for j in range(nbatches):
		if i!=j:
			rho_true[k] = np.sum(mass*rho[:,j])
			k+=1
	return L1, np.sum(rho_true)/(nbatches-1)

def Likelihood(sn):
	L = np.zeros(nbatches)
	for i in range(nbatches):
		L1, rho_true = sum_rho_true(i, sn)
		L2 = -2 * rho_true
		L[i] = L1 + L2
	return np.mean(L)


if __name__ == "__main__":
	densities_file = "./data_KL_SN/rho_estimates/MW_100M_b1_dm_part_1e6_300_KL_analysis_batch"
	file_name = "Likelihood_CM_sn.txt"

	npart = 1000
	mass = 1 / npart
	nbatches = 2
	sn_range = np.arange(0, 8, 0.2)
	num_cores= 4


	output = Parallel(n_jobs=num_cores)(delayed(Likelihood)(sn) for sn in range(len(sn_range)))
	print(output)
	write_output(file_name, output)

