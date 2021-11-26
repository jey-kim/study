	PROGRAM INPUT_SPLINE_ROTATIONS

	integer	rot_num,lat_pole,lon_pole,i,num,xstep,ystep,
     c          total_knotpoints,inside_knotpoints,boundary_knotpoints,
     c          a,b,c

	real d,e,f,g,v,w,s,t,u

	double precision rot

	open(10,file='rot.dat')
	read(10,*)

	open(11,file='boundary_spline_fit.dat')
	open(12,file='spline_fit.dat')

	read(11,*)xstep,ystep,total_knotpoints
	write(12,*)xstep,ystep,total_knotpoints


	read(11,*)inside_knotpoints
	write(12,*)inside_knotpoints

        boundary_knotpoints = total_knotpoints - inside_knotpoints      

	do i=1,boundary_knotpoints
		read(10,*)rot_num,lat_pole,lon_pole,rot
		write(12,*)rot_num,lat_pole,lon_pole,rot
	enddo
	
	close(10)
	
	read(11,*)num
	write(12,*)num
c      num = total number of midpoints		

	do i=1,num
		read(11,*)a,b,c
		read(11,*)s,t,u
		read(11,*)d,e,f,g,v,w
	
		write(12,*)
		write(12,*)a,b,c
		write(12,*)s,t,u
		write(12,*)d,e,f,g,v,w
	enddo

22	format(1x,i5,1x,i5,1x,i5,1x,i5,1x,e15.8,1x,e15.8,1x,e15.8)
23	format(1x,e10.3,1x,e10.3,1x,e10.3)
24	format(1x,e10.3,1x,e10.3,1x,e10.3,1x,f10.7,1x,f10.7,1x,f10.7)

	close(11)
	close(12)
	end
