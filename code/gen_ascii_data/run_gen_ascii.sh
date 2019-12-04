in_path="../../../MW_anisotropy/code/test_snaps/"
out_path="./"
snap_name="MWLMC3_100M_new_b1_"
nmax=20
lmax=20
rs=40.85
n_halo_part=100000000

python3 armadillo_find_bound_particles.py $in_path $snap_name "test_file" $nmax $lmax $rs $n_halo_part  --ncores=4

