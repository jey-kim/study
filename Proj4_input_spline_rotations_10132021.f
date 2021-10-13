	PROGRAM INPUT_SPLINE_ROTATIONS

	integer	a,b,c,i,num

	real e,f,g,v,w,s,t,u

	double precision d

	open(10,file='rot.dat')
	read(10,*)

	open(11,file='boundary_spline_fit.dat')
	open(12,file='spline_fit.dat')

	read(11,*)a,b,c
	write(12,*)a,b,c

c       a,b,c = X coordinate,Y coordinate, total number of rotations

	read(11,*)a
	write(12,*)a

c       a = total number of rotations - total number of boundary points

	do i=1,48
		read(10,*)a,b,c,d
		write(12,*)a,b,c,d
	enddo
	
	close(10)
	
	read(11,*)num
	write(12,*)num

c      num = total number of areas		

	do i=1,num
c		read(11,*)
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
