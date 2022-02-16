#!/bin/bash

target_location="NE"
target_date="02082022"
target_inv="L2"

target_directory="/Users/jkim/main/Joint_GNSS_INSAR/Interseismic_inversion_SCalifornia_"$target_location"_"$target_date""
target_directory_true="/Users/jkim/main/Joint_GNSS_INSAR/Interseismic_syntheticDATA_generator_SCalifornia_"$target_location"_"$target_date""

cwd=$(pwd)

echo "\n Copying files associated with the "$target_inv" inversion results from"
echo "'"$target_directory"' directory.\n" 

echo "Copying files associated with the 3-D true fields from"
echo "'"$target_directory_true"' directory.\n"
 
## NOTE !!
## for synthetic data sets: 
## The original local frame of reference for an InSAR synthetic data set
## is equal to the converted local frame.
## Thus, the true 3-D field should be recovered by the joint inversion 
## using the InSAR synthetic data in the converted frame.
 
target_file1_true="vel_vertical_scaled.gmt" # RESAMPLE THIS!
target_file2_true="vel_horizontal_cont.gmt" # RESAMPLE THIS!

cp "$target_directory_true"/"$target_file1_true" .
cp "$target_directory_true"/"$target_file2_true" .

target_file_inversion_Gmat="G_matrix.out"
target_file_inversion_GNSS_data="GNSS_input.gmt"
target_file_inversion_BC_data="BC_input.gmt"
target_file_inversion_INSAR_datasets="_data_common_new_ref.dat"

cp "$target_directory"/"$target_file_inversion_Gmat" .
cp "$target_directory"/"$target_file_inversion_GNSS_data" .
cp "$target_directory"/"$target_file_inversion_BC_data" .
cp "$target_directory"/*"$target_file_inversion_INSAR_datasets" .

DSD_InSAR_file=`ls DT*`
cp "$DSD_InSAR_file" dsd_InSAR_input.dat
rm "$DSD_InSAR_file"

ASD_InSAR_file=`ls AT*`
cp "$ASD_InSAR_file" asd_InSAR_input.dat 
rm "$ASD_InSAR_file"


for weighting_InSAR in 0.01 0.1 1 2 3 4 5 10; do
    target_file1=""$target_inv"_stat_"$weighting_InSAR""
    cp "$target_directory"/"$target_file1"* .
    target_file2=""$target_inv"_model_coef_"$weighting_InSAR""
    cp "$target_directory"/"$target_file2"* .

    target_file3=""$target_inv"_vel_horizontal_cont_pred_"$weighting_InSAR""
    cp "$target_directory"/"$target_file3"* .
    target_file4=""$target_inv"_vel_vertical_cont_pred_"$weighting_InSAR""
    cp "$target_directory"/"$target_file4"* .
done


# L2_vel_GNSS_pred_10_5_500_10.gmt
# L2_dsd_InSAR_pred_10_5_500_10.dat
# L2_asd_InSAR_pred_10_5_500_10.dat
# L2_vel_BC_pred_10_5_500_10.dat
