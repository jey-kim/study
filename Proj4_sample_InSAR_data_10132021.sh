#!/bin/bash

#STEP 2

#STEP2-a 

# Convert ipynb to py
##jupyter nbconvert Proj4_read_mat_files_10132021.ipynb --to python
python Proj4_read_mat_files_10132021.py

#STEP2-b
#interpolate data & plot data
sh Proj4_process_data_interpolation_10132021.sh


sh image_dLOS_check_Sentinel
open dLOS_masked_Sentinel_obs_1.pdf
open dLOS_masked_Sentinel_obs_2.pdf
open dLOS_masked_Sentinel_obs_3.pdf
open dLOS_masked_Sentinel_obs_4.pdf

open dLOS_interp_Sentinel_obs_1.pdf
open dLOS_interp_Sentinel_obs_2.pdf
open dLOS_interp_Sentinel_obs_3.pdf
open dLOS_interp_Sentinel_obs_4.pdf


#STEP2-c
#save coordinates of the sampled inSAR data
awk '{print $1,$2}' interpolated_A_4.dat > interpolated_coordinate.dat
awk '{print $1,$2}' A_4.dat > masked_coordinate.dat

#STEP2-d 
# Change the original frame of reference for the InSAR data to
# the same frame of reference for the basis-function responses

python Proj4_change_InSAR_ref_10132021.py


sh image_dLOS_check_Sentinel_new_ref
open dLOS_interp_Sentinel_obs_new_ref_1.pdf
open dLOS_interp_Sentinel_obs_new_ref_2.pdf
open dLOS_interp_Sentinel_obs_new_ref_3.pdf
open dLOS_interp_Sentinel_obs_new_ref_4.pdf

