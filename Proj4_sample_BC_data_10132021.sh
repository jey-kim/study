#!/bin/bash

# This code uses "spline_fit_ref.out" file.
# Velocity reference point at (-116.528 E ,33.908 N).

####################################################
## STEP 1 : getting latlong_gps.dat               ##
## Inputs : geometry.dat, sparse.log, GPS_raw.dat ## 
## Output : latlong_gps.dat                       ##
####################################################

cp geometry_hard_wiring.dat geometry.dat
# This geometry covers the entire S.California and Western Nevada
# Kim et al. 2021a & 2021b

make_sparse_geometry 
# For sparse.log file

python Proj4_make_fake_GPS_raw_for_BC_10132021.py
cp GPS_raw_fake_BC.dat GPS_raw.dat
# Make a fake GPS_raw.dat -> boundary points + reference point
# Reference point on the first row


prepare_gps_sparse2<<!
7237
!



##############################################################
## STEP 2 : obtain reference velocity value (lat, lon, rot) ##
## INPUTS : spline_fit_ref.out, latlong_gps.dat             ##
## OUTPUT : velocity.out                                    ##
##############################################################

cp spline_fit_ref.out spline_fit.out

cp latlong_gps.dat output.dat

make_x_for_output<<!
0 0 0
!

sparse_velocity<<!
Y
N
0 0 0
!

velocity 


awk '{print $1-360,$2,$3,$4,$5,$6,$7}' vel.gmt > vel_PA_ref.gmt
#cp vel.gmt vel_PA_ref.gmt

######Rotation info
awk '{ if (NR==3) print $7,$8,$9}' velocity.out > lat_long_rotation.out
#check the lat lon rot value at the reference point

#gfortran sparse_velocity_rotated.f -o sparse_velocity_rotated -ffixed-line-length-none

./sparse_velocity_rotated<<!
Y
N
!

velocity

#cp vel.gmt vel_InSAR_ref.gmt
awk '{print $1-360,$2,$3,$4,$5,$6,$7}' vel.gmt > vel_InSAR_ref.gmt
#remove the reference point from vel_UCERF.gmt

awk 'NR>1' vel_InSAR_ref.gmt > vel_InSAR_ref_no_refpoint.gmt

#plot
sh image_BC_data_check
open boundary_data_check.pdf
