		program average_strain_RECTANGULAR 
		
 			

        	open(15,file='average_strain_4_mo.out',status='old')
		
		open(77,file='average_strain_4_mo_RECTANGULAR.out')

		read(15,*)
		read(15,*)

		do i = 1, 7185
	
		read(15,*)num,rlat,rlon,exx,eyy,exy,e11,e22,e12
		read(15,*)sig_xx,sig_yy,sig_xy,se11,se22,se12
		read(15,*)tcor1,tcor2,tcor3,ycor1,ycor2,ycor3
		
		write(77,99)num,rlat,rlon,e11,e22,e12,se11,se22,se12
99    format(1x,I5,1x,2f14.5,1x,6f14.5)
                end do

		stop
		end
