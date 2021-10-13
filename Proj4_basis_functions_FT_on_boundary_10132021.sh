#!/bin/bash


make_sparse_geometry

make_raw_spline_fit_dat

cp spline_fit.dat spline_fit_raw.dat

make_x_for_fit<<!
0 0 0
!

#for errors
make_null_input_with_variance<<!
1e-8
0.25
!

cp spline_fit.dat spline_fit_null.dat

#for values  #smaller areas 
make_null_input_with_variance<<!
0.0246
0.25
!

cp spline_fit.dat spline_fit_value.dat


sed 2d GPS_raw_fake_BC.dat > GPS_raw.dat

prepare_gps_sparse2<<!
165
!

cp latlong_gps.dat latlong_FT_boundary.dat



for i in {1..140}  # Grid cells

do

for j in {1..4}    # Potentials

do


#i=5
#j=4

echo $i $j

#gfortran make_spline_fit_for_basis_fuctions.f
gfortran Proj4_make_spline_fit_for_basis_fuctions_FT.f
./a.out<<!
$i $j
!

cp spline_fit_for_basis.dat spline_fit_for_basis_"$i"_"$j".dat

cp spline_fit_for_basis.dat spline_fit.dat

make_x_for_fit<<!
0 0 0
!

sparse_fit<<!
N
Y
Y
!

cp latlong_FT_boundary.dat output.dat

make_x_for_output<<!
0 0 0
!

sparse_velocity<<!
Y
N
0 0 0
!

sparse_strain<<!
Y
N
!

velocity

cp vel.gmt vel_hori_FT_"$i"_"$j"_on_boundary.gmt
awk '{print $1-360,$2,$3,$4,$5,$6,$7}' vel.gmt > vel_hori.gmt

cp latlong_FT_boundary.dat nout.dat
strain_cont_asia_1e9
cp sigma_c.dat sigma_c_"$i"_"$j".dat
awk '{print $1,$2,$3*6.371/2}' sigma_c.dat > vel_vert.gmt
#Exx+Eyy/2 = Exx = Eyy = Ur/r
# r = the radius of the Earth
cp vel_vert.gmt vel_vert_FT_"$i"_"$j"_on_boundary.gmt


#PLOT


if [ "$j" -lt "4" ]; then

./image_basis_function_FT_hori
cp basis_function_FT_hori.ps basis_function_FT_hori_"$i"_"$j".ps
cp basis_function_FT_hori_"$i"_"$j".ps plots_FT_on_boundary
rm basis_function_FT_hori_"$i"_"$j".ps

else

./image_basis_function_FT_hori_from_vert
cp basis_function_FT_hori.ps basis_function_FT_hori_"$i"_"$j".ps
cp basis_function_FT_hori_"$i"_"$j".ps plots_FT_on_boundary
rm basis_function_FT_hori_"$i"_"$j".ps


./image_basis_function_FT_vert
cp basis_function_FT_vert.ps basis_function_FT_vert_"$i"_"$j".ps
cp basis_function_FT_vert_"$i"_"$j".ps plots_FT_on_boundary
rm basis_function_FT_vert_"$i"_"$j".ps


fi 

done
done
