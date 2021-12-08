	program make_grid

        implicit real*8 (a-h,o-z)
        open(11, file = 'map_config.txt')
        
        read(11,*)rlon_min
        read(11,*)rlon_max
        read(11,*)rlat_min
        read(11,*)rlat_max 
        read(11,*)step
        
        float_nx=(rlon_max-rlon_min)/step
        nx = nint(float_nx)
        
        float_ny=(rlat_max-rlat_min)/step
        ny = nint(float_ny)

	open(22, file = 'geometry_kin.dat')

	write(22,*)nx,ny,3

        rlong = rlon_min
	rlat = rlat_min
        ix = 0
	iy = 0

                do j = 1, nx+1 

                     write(22,*)ix,iy,1,3,3
                     write(22,25)rlat,rlong
25                   format(2f12.5)
                        rlong = rlong + step
		     ix = ix + 1

                end do 
        iy = 1
        rlat = rlat_min + step
	do i = 1,ny-1 

	ix = 0
	rlong = rlon_min 
           write(22,*)ix,iy,1,3,3
           write(22,25)rlat,rlong
           ix = ix + 1
           rlong = rlong + step
		do j = 2, nx 

		     write(22,*)ix,iy,1,3,3
		     write(22,25)rlat,rlong
		     ix = ix + 1
                        rlong = rlong + step
		
		end do
           write(22,*)ix,iy,1,3,3
           write(22,25)rlat,rlong
           ix = ix + 1

                
	    rlat = rlat + step
	    iy = iy + 1

	end do
        ix = 0
        rlong = rlon_min 
                do j = 1, nx + 1

                     write(22,*)ix,iy,1,3,3
                     write(22,25)rlat,rlong
                        rlong = rlong + step
                     ix = ix + 1
                end do

	stop

	end


