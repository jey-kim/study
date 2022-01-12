#!/bin/bash

#This is a master code


##########################################################
#    STEP 1 :   HORIZONTAL FIELD                         #
##########################################################

# cp files
#spline_fit_ref.out
#geometry_hard_wiring.dat


#make map_config.txt
#make ref_GNSS.dat


cp geometry_hard_wiring.dat geometry.dat
make_sparse_geometry

cp spline_fit_ref.out spline_fit.out

prepare_knotpoints


lon_min=$(awk '{if(NR==1) print $1}' map_config.txt)
lon_max=$(awk '{if(NR==2) print $1}' map_config.txt)
lat_min=$(awk '{if(NR==3) print $1}' map_config.txt)
lat_max=$(awk '{if(NR==4) print $1}' map_config.txt)

prepare_regular_lat_long<<!
0 0 0
2
$lon_min $lon_max $lat_min $lat_max
100 100
!


cp regular_lat_long.dat output.dat

#cp knotpoints.dat output.dat
make_x_for_output<<!
0 0 0
!

sparse_velocity<<!
Y
N
0 0 0
!

velocity
awk '{print $1-360,$2,$3,$4,$5,$6,$7}' vel.gmt > vel_UCERF3_PA.gmt
./image_UCERF_horizontal
#open UCERF3.ps

lon_ref=$(awk '{if(NR==1) print $1}' ref_GNSS.dat)
lat_ref=$(awk '{if(NR==1) print $2}' ref_GNSS.dat)
awk -v lon_ref=$lon_ref -v lat_ref=$lat_ref '{ if ($3==lon_ref && $2==lat_ref) print $7,$8,$9}' velocity.out > lat_long_rotation.out


gfortran sparse_velocity_rotated.f -o sparse_velocity_rotated -ffixed-line-length-none
./sparse_velocity_rotated<<!
Y
N
!

velocity

awk '{print $1-360,$2,$3,$4,$5,$6,$7}' vel.gmt > vel_UCERF3_corner.gmt

./image_UCERF_horizontal_corner
#open UCERF3_corner.ps

# FINAL horizontal continuous field : vel_UCERF3_corner.gmt #




##########################################################
#    STEP 2 :   VERTICAL FIELD                           #
##########################################################

# Rigid Boundary geometry
gfortran make_geometry_rigid_map_config.f
./a.out 
cp geometry_rigid_boundary.dat geometry_synthetic.dat
cp geometry_synthetic.dat geometry.dat

make_sparse_geometry

make_raw_spline_fit_dat

cp spline_fit.dat spline_fit_makeraw.dat

make_x_for_fit<<!
0 0 0
!

make_null_input_with_variance<<!
1e-8
0.25
!

cp spline_fit.dat spline_fit_null.dat


gfortran Proj4_make_spline_vertical.f
./a.out<<!
115
8
!

cp spline_fit_vertical.dat spline_fit.dat

make_x_for_fit<<!
0 0 0
!

sparse_fit<<!
Y
Y
N
!

prepare_regular_lat_long<<!
0 0 0
2
$lon_min $lon_max $lat_min $lat_max
100 100
!
cp regular_lat_long.dat output.dat

make_x_for_output<<!
0 0 0
!

cp regular_lat_long.dat nout.dat
sparse_strain<<!
Y
!

sparse_velocity<<!
Y
0 0 0
!

velocity

cp vel.gmt vel_hori_vert.gmt
awk '{print $1-360,$2,$3,$4,$5,$6,$7}' vel_hori_vert.gmt > vel_horizontal_from_vertical_displacement.gmt
strain_cont_asia_1e9

awk '{print $1,$2,$3*6.371/2}' sigma_c.dat > vel_vertical.gmt

sh image_vertical.sh
open vertical.pdf
sh image_hori_from_vert
open horizontal_from_vert.pdf



##########################################################
#    STEP 3 :   GENERATE Synthetic InSAR & GNSS data     #
##########################################################

# ARIA-tools was used to obtain the following files
# azimuthAngle_*
# incidenceAngle_*
# mask_*

track1="DT71"
gmt grd2xyz incidenceAngle_"$track1".grd > incidenceAngle_"$track1".xyz
gmt grd2xyz azimuthAngle_"$track1".grd > azimuthAngle_"$track1".xyz
gmt grd2xyz mask_"$track1".grd > mask_"$track1".xyz
python get_pointing_vector.py "$track1"

track2="AT64"
gmt grd2xyz incidenceAngle_"$track2".grd > incidenceAngle_"$track2".xyz
gmt grd2xyz azimuthAngle_"$track2".grd > azimuthAngle_"$track2".xyz
gmt grd2xyz mask_"$track2".grd > mask_"$track2".xyz
python get_pointing_vector.py "$track2"

############
############
############ catch up from here 01-12-2022
############

sh image_Px_DT.sh
sh image_Py_DT.sh
sh image_Pz_DT.sh
sh image_Px_AT.sh
sh image_Py_AT.sh
sh image_Pz_AT.sh

# resample Pi (asc and dsc) on regular lat lon
sh get_insar_data_regular_lat_lon.sh
 

python Proj4_generate_synthetic_data_12092021.py
#The vertical velocity and associated horizontal field 
# are SCALED by 0.3 !

#Errors for GNSS 5% InSAR 30 %


sh image_DT173_synthetic.sh
sh image_AT64_synthetic.sh
sh image_GNSS_synthetic.sh
open GNSS_synthetic.pdf
open dLOS_DT173_synthetic.pdf
open dLOS_AT64_synthetic.pdf
