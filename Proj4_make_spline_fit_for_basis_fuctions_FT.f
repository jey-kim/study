	program make_spline_fit_for_basis_fuctions
		
        implicit real*8 (a-h,o-z)
        
	open(10,file='spline_fit_for_basis.dat')
        open(20,file='spline_fit_raw.dat') 
	open(30,file='spline_fit_null.dat')	
	open(77,file='spline_fit_value.dat')

        read(30,*)
        read(30,*)
        read(30,*)
        read(30,*) 

        read(20,*)i,j,k
        read(20,*)l
        read(20,*)m,r1,r2,r3
        read(20,*)i1

	
	read(77,*)
	read(77,*)
	read(77,*)
	read(77,*)
	
        
	write(6,*)'enter the area number and enter the strain index'
	read(5,*)marea_index,mstrain_index
	
	write(10,*)i,j,k
	write(10,*)l
	write(10,*)m,r1,r2,r3
	write(10,*)i1

		
	do i = 1, i1
	read(30,*)	
	read(30,*)
	read(30,*)exx1,eyy1,exy1
	read(30,*)sem_xx,sem_yy,sem_xy,cc1,cc2,cc3

 	read(20,*)
	read(20,*)ii,jj,kk
	read(20,*)
	read(20,*)
	
	read(77,*)
	read(77,*)
	read(77,*)
	read(77,*)the_value_xx,the_value_yy,the_value_xy
		
	write(10,*)
	write(10,*)ii,jj,kk
		if (ii .EQ. marea_index) then 
			if (mstrain_index .EQ. 1) then
				write(10,23)the_value_xx,eyy1,exy1
			        write(10,24)sem_xx,sem_yy,sem_xy,cc1,cc2,cc3
			else if (mstrain_index .EQ. 2) then
				write(10,23)exx1,the_value_yy,exy1
			        write(10,24)sem_xx,sem_yy,sem_xy,cc1,cc2,cc3
			else if (mstrain_index .EQ. 3) then
				write(10,23)exx1,eyy1,the_value_xy
				write(10,24)sem_xx,sem_yy,sem_xy,cc1,cc2,cc3
			else if (mstrain_index .EQ. 4) then
				write(10,23)the_value_xx/6.371,the_value_yy/6.371,exy1
				write(10,24)sem_xx,sem_yy,sem_xy,cc1,cc2,cc3
			end if
		else
			write(10,23)exx1,eyy1,exy1
			write(10,24)sem_xx,sem_yy,sem_xy,cc1,cc2,cc3
		end if

23	format(3e15.5)
24	format(6e15.5)
	
	end do

		
	stop
	end
