#!/bin/bash

# Step 1: Copy InSAR data to the work directory
sh Proj4_cp_inSAR_BC_data_10132021.sh


# Step 2: Resample InSAR data and perform interpolation ?
sh Proj4_sample_InSAR_data_10132021.sh


# Step 3: Resample UCERF3 boundary condition velocity data.
#         The insar data has a reference point at (-116.528; 33.908; coherent pixel as a ref point)
sh Proj4_sample_BC_data_10132021.sh


# Step 4: Obtain basis functions : 
#         Force Balance equations inside & boundary basis function (from 3 rotations :x,y,and z)
sh Proj4_generate_basis_functions_10132021.sh 


# Step 5: Invert them. 
#         Weighting factor for InSAR will be changed from 1 to 8.
#         An weighting factor will be inversely taken by the algorithm (e.g., 8 => 1/8)
sh Proj4_inversion_10132021.sh 

