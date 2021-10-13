#!/bin/bash

cp D_1.dat A_3.dat #Descending 1
cp D_2.dat A_4.dat #Descending 2

for i in {1..4}; do

awk '{print $1,$2,$6}' A_"$i".dat > input.gmt
sh image_dLOS_interpolation_input_Sentinel
cp input.xyz rate_A_"$i".dat
rm input.gmt input.xyz

awk '{print $1,$2,$3}' A_"$i".dat > input.gmt
sh image_dLOS_interpolation_input_Sentinel
cp input.xyz Px_A_"$i".dat
rm input.gmt input.xyz

awk '{print $1,$2,$4}' A_"$i".dat > input.gmt
sh image_dLOS_interpolation_input_Sentinel
cp input.xyz Py_A_"$i".dat
rm input.gmt input.xyz

awk '{print $1,$2,$5}' A_"$i".dat > input.gmt
sh image_dLOS_interpolation_input_Sentinel
cp input.xyz Pz_A_"$i".dat
rm input.gmt input.xyz

paste rate_A_"$i".dat Px_A_"$i".dat Py_A_"$i".dat Pz_A_"$i".dat > final.dat
awk '{print $1,$2,$3,$6,$9,$12}' final.dat > interpolated_A_"$i".dat
rm final.dat

done


