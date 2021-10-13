#!/bin/bash 

rm stat.txt

weighting_BC=1
for weighting_InSAR in {1..20}; do
 
python Proj4_inversion_BC_and_InSAR_weighting_10132021.py "$weighting_BC" "$weighting_InSAR" >> stat.txt

sh image_BC_predicted_final
cp boundary_data_and_prediction_final.pdf boundary_data_and_prediction_final_"$weighting_InSAR"_"$weighting_BC".pdf 

sh image_dLOS_predicted_final
cp dLOS_pred_A_1_final.pdf dLOS_pred_A_1_final_"$weighting_InSAR"_"$weighting_BC".pdf
cp dLOS_pred_A_2_final.pdf dLOS_pred_A_2_final_"$weighting_InSAR"_"$weighting_BC".pdf
cp dLOS_pred_A_3_final.pdf dLOS_pred_A_3_final_"$weighting_InSAR"_"$weighting_BC".pdf
cp dLOS_pred_A_4_final.pdf dLOS_pred_A_4_final_"$weighting_InSAR"_"$weighting_BC".pdf

done
