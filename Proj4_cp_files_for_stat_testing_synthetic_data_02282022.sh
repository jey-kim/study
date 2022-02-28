#!/bin/bash

## 1. ## Define Params.

cp_sysnthetic_field="True"
target_path="/Users/jkim/main/Joint_GNSS_INSAR/Interseismic_inversion_SCalifornia"
target_path_syn="/Users/jkim/main/Joint_GNSS_INSAR/Interseismic_syntheticDATA_generator_SCalifornia"
target_loc="SE"
target_date="02172022"
target_date_syn="02082022"
target_inv="L2"
loop_param_dir="first_batch second_batch third_batch"
weighting_params="1.406_0.156" #wI_wG


## 2. ## Define more Params. 
cp_file_Gmat="G_matrix.out"
cp_file_GNSS_data="GNSS_input.gmt"
cp_file_BC_data="BC_input.gmt"
cp_file_INSAR_data="_data_common_new_ref.dat"
cp_file_syn_vrt="vel_vertical_scaled.gmt" 
cp_file_syn_hrz="vel_horizontal_cont.gmt" 


## 3. ## Create a subdir where all the inversion results will be saved
mkdir -p input_files_"$target_inv"

## 4. ## Save variables for directiory names  
cp_dir=""$target_path"_"$target_loc"_"$target_date"_"$target_inv""
cp_dir_syn=""$target_path_syn"_"$target_loc"_"$target_date_syn""
cwd=$(pwd)
save_dir=""$cwd"/input_files_"$target_inv""

echo "\n Copying files associated with the "$target_inv" inversion results from"
echo " '"$cp_dir"_ith_batch' directory,"
echo "   where i in "$loop_param_dir".\n" 
echo "Copying files associated with the synthetic true fields from"
echo " '"$cp_dir_syn"' directory.\n"
 

## 5. ## CP synthetic surface 3-D field
## !! NOTE !! for synthetic data sets: 
## The original local frame of reference for an InSAR synthetic data set
## is equal to the converted local frame.
## Thus, the true 3-D field should be recovered by the joint inversion 
## using the InSAR synthetic data in the converted frame.
if [[ $cp_sysnthetic_field == "True" ]]; then
    cp "$cp_dir_syn"/"$cp_file_syn_vrt" "$save_dir"
    cp "$cp_dir_syn"/"$cp_file_syn_hrz" "$save_dir"
fi

## 6. ## CP data files and Gmatrix file
fir="$(echo $loop_param_dir | awk '{print $1}')"

cp "$cp_dir"_"$fir"/"$cp_file_Gmat" "$save_dir"
cp "$cp_dir"_"$fir"/"$cp_file_GNSS_data" "$save_dir"
cp "$cp_dir"_"$fir"/"$cp_file_BC_data" "$save_dir"
cp "$cp_dir"_"$fir"/*"$cp_file_INSAR_data" "$save_dir"

cd "$save_dir"
DSD_InSAR_file=`ls DT*`
cp "$DSD_InSAR_file" dsd_InSAR_input.dat
rm "$DSD_InSAR_file"
ASD_InSAR_file=`ls AT*`
cp "$ASD_InSAR_file" asd_InSAR_input.dat 
rm "$ASD_InSAR_file"
cd ../

## 7. ## CP inversion results

for batch in $loop_param_dir; do 
        
    cp_file1=""$target_inv"_stat_"$weighting_params""
    cp_file2=""$target_inv"_model_coef_"$weighting_params""
    cp_file3=""$target_inv"_vel_horizontal_cont_pred_"$weighting_params""
    cp_file4=""$target_inv"_vel_vertical_cont_pred_"$weighting_params""
    
    cp "$cp_dir"_"$batch"/"$cp_file1"* "$save_dir"
    cp "$cp_dir"_"$batch"/"$cp_file2"* "$save_dir"
    cp "$cp_dir"_"$batch"/"$cp_file3"* "$save_dir"
    cp "$cp_dir"_"$batch"/"$cp_file4"* "$save_dir"        

done


echo "DONE!"
