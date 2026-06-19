import os

import numpy as np


"""
#Non symmetrized, no covariance
1. write header properly (for now assume some value for nbins, but will be generated in analysis script)
2. Bring numerical data in ALF format
3. tau 0 data value should be some value with some large error, so that relative error becomes large? or something like this idk. should also be non zero value.
3. run and save as g_dat file.


-manually-
4. cp parameters file from script analysis_start
5. Overwrite parameter file with correct parameters (do manually for now by taking from parameter file)
6. Test run SAC -> did run

7. Compare with result from Superlattice SAC

#Non symmetrized, with covariance
1. now use the same file from before and write the covariance matrix row by row in the g_dat file
2. change ncov parameter
3. Run SAC

#In analysis script:
1. Implement symmetrization before calculating the covariance matrix, make sure to symmetrize correctly, ignoring the 0 value and so on..
2. Now use the symmetrized data with the alf port
3. now run again.

#NOW verify channel?
"""



#TODO: Number of bins should not be fixed
def generate_header_string(param_path: str, nbins: int = 100000) -> str:
    norb = 1
    channel = "PH_C"

    def return_beta_nt(param_path: str):
        with open(param_path, "r") as file_handle:
            lines = [line.strip() for line in file_handle if line.strip()]

        beta = lines[lines.index("beta") + 1]
        ntau = lines[lines.index("Nt_intervals_full") + 1]
        return beta, ntau
    
    beta, ntau = return_beta_nt(param_path)
    
    dtau=float(beta)/float(ntau)

    return dtau, f"{int(ntau):>12}{nbins:>11}  {float(beta):.17E}{norb:>11} {channel}\n"

def load_data(dtau, data_path: str):
    data = np.loadtxt(data_path)
    data[0,2] = 1e10
    data[:,0]=data[:,0]*dtau
    return data

def load_cov(cov_path: str):
    return np.loadtxt(cov_path)


def generate_g_dat_file(save_path: str, header: str, data: np.array, cov: np.array):
    g_dat_path = os.path.join(save_path, "g_dat")
    
    with open(g_dat_path, "w") as file_handle:
        file_handle.write(header)
        if not header.endswith("\n"):
            file_handle.write("\n")
        np.savetxt(file_handle, np.asarray(data), fmt="%.17E")
        np.savetxt(file_handle, np.asarray(cov).ravel(order="C"), fmt="%.17E")



if __name__=="__main__":
    base_dir="U_5.5_BETA_2.0_NT_256"
    test_file=os.path.join(base_dir,"K1_0_K2_0_S1_0_S2_0_L1_0_L2_0.txt")
    PARAM_PATH=os.path.join(base_dir,"parameter_an_cont_ED_corr0.txt")
    COV_PATH=os.path.join(base_dir,"K1_0_K2_0_S1_0_S2_0_L1_0_L2_0_covariance_real.txt")
    save_path=os.path.join(base_dir,"TEST_PORT")
    ALF_DIR="/Users/luis/Documents/_WORK/Coding/ALF_CLEAN/ALF"
    
    dtau, header=generate_header_string(PARAM_PATH)
    data=load_data(dtau, test_file)
    cov=load_cov(COV_PATH)
    generate_g_dat_file(save_path, header, data, cov)
