#!/bin/bash

# This code generates basis functions for BC and InSAR data points
# See the equation for the G matrix elements! 
# (This equation was written in 9/24/2021)
# Jeonghyeop.kim@gmail.com
# 09/24/2021

######################################################################
#                                                                   ##
# STEP 4-a : Prepare latlong_gps.dat for the boundary points        ##
#          : as well as for the InSAR data points.                  ##
#  OUTPUTS : latlong_BC.dat & latlong_InSAR.dat                     ##
#                                                                   ##
#                                                                   ##
#                                                                   ##
######################################################################


cp geometry_BC.dat geometry.dat
make_sparse_geometry
make_raw_spline_fit_dat
make_x_for_fit<<!
0 0 0
!

### STEP 4-a-1 : latlong_gps.dat for BOUNDARY POINTS
##Delete the second line (1st line: header; 2nd line: InSAR ref point)
sed 2d GPS_raw_fake_BC.dat > GPS_raw.dat

prepare_gps_sparse2<<!
165
!

cp latlong_gps.dat latlong_BC_boundary.dat


### STEP 4-a-b: latlong_gps.dat for InSAR DATA POINTS
## make fake GPS_raw.dat to generate latlong_gps.dat for InSAR data points 
python Proj4_make_fake_GPS_raw_for_InSAR_10132021.py

cp GPS_raw_fake_InSAR.dat GPS_raw.dat

prepare_gps_sparse2<<!
165
!

cp latlong_gps.dat latlong_BC_InSAR.dat


######################################################################
#STEP 4-b : Generate basis functions related to boundary conditions ##
#         : at Boundary points as well as at InSAR data points.     ##
# OUTPUTS : vel_BC_x_"$i"_on_boundary.gmt                           ##
#         : vel_BC_y_"$i"_on_boundary.gmt                           ##
#         : vel_BC_z_"$i"_on_boundary.gmt                           ##
#         : vel_BC_x_"$i"_on_InSAR.gmt                              ##
#         : vel_BC_y_"$i"_on_InSAR.gmt                              ##
#         : vel_BC_z_"$i"_on_InSAR.gmt                              ##
#                                                                   ##
# PLOTS   : see "/plots_BC_on_boundary" & /plots_BC_on_InSAR        ##
######################################################################

sh Proj4_basis_functions_BC_on_boundary_10132021.sh 
sh Proj4_basis_functions_BC_on_InSAR_10132021.sh 

######################################################################
#STEP 4-c : Generate basis functions related to force terms         ##
#         : at Boundary points as well as at InSAR data points.     ##
# OUTPUTS : vel_FT_exx_"$i"_on_boundary.gmt                         ##
#         : vel_FT_eyy_"$i"_on_boundary.gmt                         ##
#         : vel_FT_exy_"$i"_on_boundary.gmt                         ##
#         : vel_FT_vtc_"$i"_on_boundary.gmt    vertical  exx=eyy    ##
#         : vel_FT_exx_"$i"_on_InSAR.gmt                            ##
#         : vel_FT_eyy_"$i"_on_InSAR.gmt                            ##
#         : vel_FT_exy_"$i"_on_InSAR.gmt                            ##
#         : vel_FT_vtc_"$i"_on_InSAR.gmt       vertical  exx=eyy    ##
# PLOTS   : see "/plots_FT_on_boundary" & /plots_FT_on_InSAR        ##
######################################################################
######################################################################

#cp geometry_free_deforming.dat geometry.dat
#Bill suggested to use geometry.dat with rigid boundaries! 
# 10/6/2021
#gfortran make_geometry_generic_rigid_boundary_INSAR.f
#./a.out

cp geometry_rigid_boundary.dat geometry.dat


sh Proj4_basis_functions_FT_on_boundary_10132021.sh
sh Proj4_basis_functions_FT_on_InSAR_10132021.sh
