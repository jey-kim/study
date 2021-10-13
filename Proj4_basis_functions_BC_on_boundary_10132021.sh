#!/bin/bash


cp latlong_BC_boundary.dat output.dat

make_x_for_output<<!
0 0 0
!

for i in `seq -f "%03g" 1 24`; do


#x-rot
cp xrot"$i".dat rot.dat

#gfortran input_spline_rotations.f
gfortran Proj4_input_spline_rotations_10132021.f
./a.out

sparse_fit <<!
N
Y
N
!

sparse_velocity<<!
Y
0 0 0
!
velocity
cp vel.gmt vel_BC_x_"$i"_on_boundary.gmt
awk '{print $1-360,$2,$3,$4,$5,$6,$7}' vel.gmt > vel_input.gmt

./image_basis_function_BC
cp basis_function_BC.ps basis_function_BC_x_"$i"_on_boundary.ps
cp basis_function_BC_x_"$i"_on_boundary.ps plots_BC_on_boundary
rm basis_function_BC_x_"$i"_on_boundary.ps


#y-rot

cp yrot"$i".dat rot.dat

#gfortran input_spline_rotations.f
gfortran Proj4_input_spline_rotations_10132021.f
./a.out

sparse_fit <<!
N
Y
N
!

sparse_velocity<<!
Y
0 0 0
!
velocity
cp vel.gmt vel_BC_y_"$i"_on_boundary.gmt
awk '{print $1-360,$2,$3,$4,$5,$6,$7}' vel.gmt > vel_input.gmt
./image_basis_function_BC
cp basis_function_BC.ps basis_function_BC_y_"$i"_on_boundary.ps
cp basis_function_BC_y_"$i"_on_boundary.ps plots_BC_on_boundary
rm basis_function_BC_y_"$i"_on_boundary.ps

#z-rot

cp zrot"$i".dat rot.dat

#gfortran input_spline_rotations.f
gfortran Proj4_input_spline_rotations_10132021.f
./a.out

sparse_fit <<!
N
Y
N
!

sparse_velocity<<!
Y
0 0 0
!
velocity
cp vel.gmt vel_BC_z_"$i"_on_boundary.gmt
awk '{print $1-360,$2,$3,$4,$5,$6,$7}' vel.gmt > vel_input.gmt
./image_basis_function_BC
cp basis_function_BC.ps basis_function_BC_z_"$i"_on_boundary.ps
cp basis_function_BC_z_"$i"_on_boundary.ps plots_BC_on_boundary
rm basis_function_BC_z_"$i"_on_boundary.ps

done


