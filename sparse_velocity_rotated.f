       
       program velstra
       implicit real*8 (a-h,o-z)
       dimension ylat(100000),ylong(100000),vary(3,3,100000),
	1	vely(3,100000)
        
	open(unit=1,file='x_for_output.dat',form='unformatted',status='old')
	read(1) ny
          call velocity(ny,ylat,ylong,vely,vary)
           open(30,file='velocity.out')
         write(30,34)
34       format(1x,'  no.  ',2x,'   LAT.   ',2x,'   LONG.  ',2x,
	1	'     Vlong     ',2x,'     Vlat      ',2x,
	2	'     Vrot      ',2x,' Pole LAT.',2x,
	3	'Pole LONG.',2x,'     Vrot      ')
           write(30,*)' '
	conv=45.0d0/datan(1.0d0)
             do i=1,ny
		alat=ylat(i)/conv
		along=ylong(i)/conv
		coslat=dcos(alat)
		sinlat=dsin(alat)
		coslong=dcos(along)
		sinlong=dsin(along)
		rot1=sinlong*vely(2,i)
	1		+coslong*(-sinlat*vely(1,i)+coslat*vely(3,i))
		rot2=-coslong*vely(2,i)
	1		+sinlong*(-sinlat*vely(1,i)+coslat*vely(3,i))
		rot3=coslat*vely(1,i)+sinlat*vely(3,i)
		prot=dsqrt(rot1**2+rot2**2+rot3**2)
		if (prot.gt.0.0d0) then
			plat=conv*dasin(rot3/prot)
		else
			plat=0.0d0
		end if
		if ((rot1.ne.0.0d0).or.(rot2.ne.0.0d0)) then
			plong=conv*datan2(rot2,rot1)
		else
			plong=0.0d0
		end if
              write(30,33)i,ylat(i),ylong(i),vely(1,i),vely(2,i),
	1	vely(3,i),plat,plong,prot
33      format(1x,i7,2x,f10.5,2x,f10.5,2x,e15.8,2x,e15.8,2x,e15.8,
	1	2x,f10.5,2x,f10.5,2x,e15.8)
             end do
         write(30,*)' '
         write(30,*)'Variance Covariance Matrices'
         write(30,*)' '
         do i=1,ny
         write(30,36)i,ylat(i),ylong(i)
          do k=1,3
             write(30,37)(vary(k,j,i),j=1,3)
          end do
          write(30,*)' '
         end do
36       format(1x,'no.',i7,2x,'LAT.=',f10.5,2x,'LONG.=',f10.5)      
37       format(1x,e15.8,2x,e15.8,2x,e15.8)
c         stop
         end
c
c
c
c
	SUBROUTINE VELOCITY(ny,rlat,rlong,vely,vary)
c
c  Input: (1) ny = number of points
c         (4) rlat(1:ny) = standard latitude in degrees of the ny points
c         (5) rlong(1:ny) = standard longitude in degrees of the ny points
c
c  Output: (6) vely(1:3,1:ny) = the 3-component velocity vectors at the ny
c             points. The first component is the velocity in the direction of
c             increasing longitude, the second is the velocity in the
c             direction of increasing latitude, and the third component is the
c             rate of rotation about the radial axis
c          (7) vary(1:3,1:3,1:ny) = the corresponding 3x3 variance-covariance
c             matrices at the ny points
c
c
c
	implicit real*8 (a-h,o-z)
	parameter (nlong_max=301,nlat_max=370,nval_max=24500,
	1	nx_max=3*nval_max,nxx_max=35000000)
c
	dimension rlat(ny),rlong(ny),vely(3,ny),vary(3,3,ny)
	dimension xcap(3),dxxcap(3),dyxcap(3),a1(3,0:1,0:1,0:1,0:1),
	1	a2(3,0:1,0:1,0:1,0:1),a3(3,0:1,0:1,0:1,0:1),
	1	ival(0:nlong_max,0:nlat_max),idel(0:nlong_max,0:nlat_max),
	2	val(3,nval_max),w0(3),
	3	a(3,48),x(nx_max),ivala(-1:2,-1:2),		! new
	7	lena(nx_max),sdinv(nxx_max),asd(nx_max,3)
	character*1 ans,yes,no
	data yes,no/'Y','N'/
	equivalence (val,x)
c
	open(unit=43,file='spline_fit.out',form='unformatted',status='old')
	read(43) nlong,nlat

	if (nlong.gt.nlong_max) 
	1	stop 'AS GREATER THAN nlong_max, ARRAYS HAVE TO BE ENLARGED'
	if (nlat.gt.nlat_max) 
	1	stop 'AS GREATER THAN nlat_max, ARRAYS HAVE TO BE ENLARGED'

	do j=0,nlat
		read(43) (ival(i,j),i=0,nlong)
	end do
	do j=0,nlat
		read(43) (idel(i,j),i=0,nlong)
	end do
	read(43) nval,varx,nrect,kvar

	if (nval.gt.nval_max) 
	1	stop 'AS GREATER THAN nval_max, ARRAYS HAVE TO BE ENLARGED'

	read(43) ((val(i,j),i=1,3),j=1,nval)
	print *,'Number of independent rotations =',nval
	print *,' Upper bound on RMS velocity =',varx
	print *,'  Number of rectangles =',nrect
	if (kvar.eq.0) print *,'   No variances and covariances'
10	print *,'Do you want to continue (Y/N)?'
	read(*,501) ans
	if (ans.eq.no) stop
	if (ans.ne.yes) goto 10
	if (kvar.eq.0) then
		ans=no
	else
20		print *,'Do you want the variance-covariance matrices (Y/N)?'
		read(*,501) ans
		if ((ans.ne.yes).and.(ans.ne.no)) goto 20
	end if
	nx=3*nval
	if (ans.eq.yes) then
		read(43) nxfree
		ixxmax=0
		do ix=1,nxfree
			read(43) lena(ix)
			ixxmin=ixxmax+1
			ixxmax=lena(ix)
			read(43) (sdinv(ixx),ixx=ixxmin,ixxmax)
		end do
		nxx=lena(nxfree)
	end if
	print *,'Enter the latidude and longitude (in degrees) of the Euler'
	print *,' pole for the frame of reference and the rotation rate to'
	print *,'  be removed'
	open(25,file='lat_long_rotation.out',status='old')
	read(25,*)plat,plong,prot
	conv=datan(1.0d0)/45.0d0
	plat=conv*plat
	plong=conv*plong
	w0(1)=prot*dcos(plat)*dcos(plong)
	w0(2)=prot*dcos(plat)*dsin(plong)
	w0(3)=prot*dsin(plat)
	do i=1,nval
	do k=1,3
		val(k,i)=val(k,i)-w0(k)
	end do
	end do

	do iy=1,ny
c             print*,iy
		read(1) xcoord,ycoord,rlong(iy),rlat(iy)
		read(1) (xcap(k),k=1,3)
		read(1) (dxxcap(k),k=1,3)
		read(1) (dyxcap(k),k=1,3)
		i0=min0(idint(xcoord),nlong-1)
		j0=min0(idint(ycoord),nlat-1)
		i1=i0+1
		j1=j0+1
		xcoord=xcoord-i0
		ycoord=ycoord-j0
		call COEFFS(xcoord,ycoord,a1,a2,a3,
	1		xcap,dxxcap,dyxcap)
		call RELATE(i0,i1,j0,j1,nlong,nlat,
	1		ival,idel,a1,a2,a3,nx,x,ivala,a,vely(1,iy),	! new
	2		nlong_max,nlat_max,nx_max)
		if (ans.eq.yes) then
c                     write(*,*) 'Calculating variances for iy=',iy
			do l=1,3
			do k=1,nx
				asd(k,l)=0.0d0
			end do
			end do
			ixmin=nxfree+1
			do ija=0,15
				ia=mod(ija,4)-1
				ja=ija/4-1
				if (ivala(ia,ja).ne.0) then
					ija1=1+3*ija
					ija2=ija1+1
					ija3=ija2+1
					i=i0+ia
					j=j0+ja
					ixa1=1+3*(ival(i,j)-1)
					ixa2=ixa1+1
					ixa3=ixa2+1
					do l=1,3
						asd(ixa1,l)=asd(ixa1,l)
	1						+a(l,ija1)
						asd(ixa2,l)=asd(ixa2,l)
	1						+a(l,ija2)
						asd(ixa3,l)=asd(ixa3,l)
	1						+a(l,ija3)
					end do
					ixmin=min0(ixmin,ixa1)
				end if
			end do
			do l=1,3
				call INV(ixmin,nxfree,nxx,lena,sdinv,asd(1,l))
			end do
			do j=1,3
			do i=1,j
				s=0.0d0
				do k=ixmin,nxfree
					s=s+asd(k,i)*asd(k,j)
				end do
				vary(i,j,iy)=s
				vary(j,i,iy)=s
			end do
			end do
		end if
	end do
	return
c
501	format(a)
	end
c
c
c
c
	SUBROUTINE RELATE(i0,i1,j0,j1,nlong,nlat,
	1	ival,idel,a1,a2,a3,nx,x,ivala,a,y,	! new
	2	nlong_max,nlat_max,nx_max)
c
	implicit real*8 (a-h,o-z)
	dimension ival(0:nlong_max,0:nlat_max),idel(0:nlong_max,0:nlat_max),
	2	a1(3,0:1,0:1,0:1,0:1),a2(3,0:1,0:1,0:1,0:1),
	3	a3(3,0:1,0:1,0:1,0:1),x(nx_max),ivala(-1:2,-1:2),a(3,48),y(3),
	4	idela(-1:2,-1:2)	! new
c
	do ix=1,48
	do iy=1,3
		a(iy,ix)=0.0d0
	end do
	end do
c
	w=1.0d0
	i0=i1-1
	j0=j1-1

	do ija=0,15
		ia=mod(ija,4)-1
		ja=ija/4-1
		ivala(ia,ja)=ija+1
	end do
	if (i0.eq.0) then
		do ja=-1,2
			ivala(-1,ja)=0
		end do
	end if
	if (i1.eq.nlong) then
		do ja=-1,2
			ivala(2,ja)=0
		end do
	end if
	if (j0.eq.0) then
		do ia=-1,2
			ivala(ia,-1)=0
		end do
	end if
	if (j1.eq.nlat) then
		do ia=-1,2
			ivala(ia,2)=0
		end do
	end if
	do ija=0,15
		ia=mod(ija,4)-1
		ja=ija/4-1
		if (ivala(ia,ja).ne.0) then
			i=i0+ia
			j=j0+ja
			if (ival(i,j).eq.0) then
				ivala(ia,ja)=0
			else
				idela(ia,ja)=idel(i,j)
			end if
		end if
	end do
	

	if (j0.eq.0) then
		if (i0.eq.0)
	1		call A00(w,0,0,ivala,idela,
	2			a1(1,0,0,0,0),a2(1,0,0,0,0),
	3			a3(1,0,0,0,0),a)
		if (i0.gt.0)
	1		call A_0(w,0,0,ivala,idela,
	2			a1(1,0,0,0,0),a2(1,0,0,0,0),
	3			a3(1,0,0,0,0),a)
		if (i1.lt.nlong)
	1		call A_0(w,1,0,ivala,idela,
	2			a1(1,0,0,1,0),a2(1,0,0,1,0),
	3			a3(1,0,0,1,0),a)
		if (i1.eq.nlong)
	1		call AN0(w,1,0,ivala,idela,
	2			a1(1,0,0,1,0),a2(1,0,0,1,0),
	3			a3(1,0,0,1,0),a)
	end if
	if (j0.gt.0) then
		if (i0.eq.0)
	1		call A0_(w,0,0,ivala,idela,
	2			a1(1,0,0,0,0),a2(1,0,0,0,0),
	3			a3(1,0,0,0,0),a)
		if (i0.gt.0)
	1		call A__(w,0,0,ivala,idela,
	2			a1(1,0,0,0,0),a2(1,0,0,0,0),
	3			a3(1,0,0,0,0),a)
		if (i1.lt.nlong)
	1		call A__(w,1,0,ivala,idela,
	2			a1(1,0,0,1,0),a2(1,0,0,1,0),
	3			a3(1,0,0,1,0),a)
		if (i1.eq.nlong)
	1		call AN_(w,1,0,ivala,idela,
	2			a1(1,0,0,1,0),a2(1,0,0,1,0),
	3			a3(1,0,0,1,0),a)
	end if
	if (j1.lt.nlat) then
		if (i0.eq.0)
	1		call A0_(w,0,1,ivala,idela,
	2			a1(1,0,0,0,1),a2(1,0,0,0,1),
	3			a3(1,0,0,0,1),a)
		if (i0.gt.0)
	1		call A__(w,0,1,ivala,idela,
	2			a1(1,0,0,0,1),a2(1,0,0,0,1),
	3			a3(1,0,0,0,1),a)
		if (i1.lt.nlong)
	1		call A__(w,1,1,ivala,idela,
	2			a1(1,0,0,1,1),a2(1,0,0,1,1),
	3			a3(1,0,0,1,1),a)
		if (i1.eq.nlong)
	1		call AN_(w,1,1,ivala,idela,
	2			a1(1,0,0,1,1),a2(1,0,0,1,1),
	3			a3(1,0,0,1,1),a)
	end if
	if (j1.eq.nlat) then
		if (i0.eq.0)
	1		call A0N(w,0,1,ivala,idela,
	2			a1(1,0,0,0,1),a2(1,0,0,0,1),
	3			a3(1,0,0,0,1),a)
		if (i0.gt.0)
	1		call A_N(w,0,1,ivala,idela,
	2			a1(1,0,0,0,1),a2(1,0,0,0,1),
	3			a3(1,0,0,0,1),a)
		if (i1.lt.nlong)
	1		call A_N(w,1,1,ivala,idela,
	2			a1(1,0,0,1,1),a2(1,0,0,1,1),
	3			a3(1,0,0,1,1),a)
		if (i1.eq.nlong)
	1		call ANN(w,1,1,ivala,idela,
	2			a1(1,0,0,1,1),a2(1,0,0,1,1),
	3			a3(1,0,0,1,1),a)
	end if

	do iy=1,3
		y(iy)=0.0d0
	end do
	do ija=0,15
		ia=mod(ija,4)-1
		ja=ija/4-1
		if (ivala(ia,ja).ne.0) then
			ija1=1+3*ija
			ija2=ija1+1
			ija3=ija2+1
			i=i0+ia
			j=j0+ja
			ixa1=1+3*(ival(i,j)-1)
			ixa2=ixa1+1
			ixa3=ixa2+1
			do iy=1,3
				y(iy)=y(iy)+a(iy,ija1)*x(ixa1)
	1				+a(iy,ija2)*x(ixa2)
	2				+a(iy,ija3)*x(ixa3)
			end do
		end if
	end do

	return
	end
c
c
c
c
	SUBROUTINE A00(w,i,j,ival,idel,a1,a2,a3,a)
c
	implicit real*8 (a-h,o-z)
	dimension ival(-1:2,-1:2),idel(-1:2,-1:2),
	1	a1(3,0:1,0:1),a2(3,0:1,0:1),a3(3,0:1,0:1),
	2	a(3,48),DX(0:1,-1:1),DY(0:1,-1:1)
	data DX/0.0,0.0,	1.0,-1.0,	0.0,1.0/
	data DY/0.0,0.0,	1.0,-1.0,	0.0,1.0/
c
	im=i
	ip=i+1
	jm=j
	jp=j+1
	call CORE(w,i,j,ival,idel,im,ip,jm,jp,DX,DY,a1,a2,a3,a)
	return
	end
c
c
c
c
	SUBROUTINE A_0(w,i,j,ival,idel,a1,a2,a3,a)
c
	implicit real*8 (a-h,o-z)
	dimension ival(-1:2,-1:2),idel(-1:2,-1:2),
	1	a1(3,0:1,0:1),a2(3,0:1,0:1),a3(3,0:1,0:1),
	2	a(3,48),DX(0:1,-1:1),DY(0:1,-1:1)
	data DX/0.0,-0.5,	1.0,0.0,	0.0,0.5/
	data DY/0.0,0.0,	1.0,-1.0,	0.0,1.0/
	im=i-1
	ip=i+1
	jm=j
	jp=j+1
	if (ival(im,j).eq.0) then
		call A00(w,i,j,ival,idel,a1,a2,a3,a)
		return
	end if
	if (ival(ip,j).eq.0) then
		call AN0(w,i,j,ival,idel,a1,a2,a3,a)
		return
	end if
	if (ival(im,jp).eq.0) then
		call A00(w,i,j,ival,idel,a1,a2,a3,a)
		return
	end if
	if (ival(ip,jp).eq.0) then
		call AN0(w,i,j,ival,idel,a1,a2,a3,a)
		return
	end if
	call CORE(w,i,j,ival,idel,im,ip,jm,jp,DX,DY,a1,a2,a3,a)
	return
	end
c
c
c
c
	SUBROUTINE AN0(w,i,j,ival,idel,a1,a2,a3,a)
c
	implicit real*8 (a-h,o-z)
	dimension ival(-1:2,-1:2),idel(-1:2,-1:2),
	1	a1(3,0:1,0:1),a2(3,0:1,0:1),a3(3,0:1,0:1),
	2	a(3,48),DX(0:1,-1:1),DY(0:1,-1:1)
	data DX/0.0,-1.0,	1.0,1.0,	0.0,0.0/
	data DY/0.0,0.0,	1.0,-1.0,	0.0,1.0/
c
	im=i-1
	ip=i
	jm=j
	jp=j+1

	call CORE(w,i,j,ival,idel,im,ip,jm,jp,DX,DY,a1,a2,a3,a)
	return
	end
c
c
c
c
	SUBROUTINE A0_(w,i,j,ival,idel,a1,a2,a3,a)
c
	implicit real*8 (a-h,o-z)
	dimension ival(-1:2,-1:2),idel(-1:2,-1:2),
	1	a1(3,0:1,0:1),a2(3,0:1,0:1),a3(3,0:1,0:1),
	2	a(3,48),DX(0:1,-1:1),DY(0:1,-1:1)
	data DX/0.0,0.0,	1.0,-1.0,	0.0,1.0/
	data DY/0.0,-0.5,	1.0,0.0,	0.0,0.5/
c
	im=i
	ip=i+1
	jm=j-1
	jp=j+1
	if (ival(i,jm).eq.0) then
		call A00(w,i,j,ival,idel,a1,a2,a3,a)
		return
	end if
	if (ival(i,jp).eq.0) then
		call A0N(w,i,j,ival,idel,a1,a2,a3,a)
		return
	end if
	if (ival(ip,jm).eq.0) then
		call A00(w,i,j,ival,idel,a1,a2,a3,a)
		return
	end if
	if (ival(ip,jp).eq.0) then
		call A0N(w,i,j,ival,idel,a1,a2,a3,a)
		return
	end if
	call CORE(w,i,j,ival,idel,im,ip,jm,jp,DX,DY,a1,a2,a3,a)
	return
	end
c
c
c
c
	SUBROUTINE A__(w,i,j,ival,idel,a1,a2,a3,a)
c
	implicit real*8 (a-h,o-z)
	dimension ival(-1:2,-1:2),idel(-1:2,-1:2),
	1	a1(3,0:1,0:1),a2(3,0:1,0:1),a3(3,0:1,0:1),
	2	a(3,48),DX(0:1,-1:1),DY(0:1,-1:1)
	data DX/0.0,-0.5,	1.0,0.0,	0.0,0.5/
	data DY/0.0,-0.5,	1.0,0.0,	0.0,0.5/
c
	im=i-1
	ip=i+1
	jm=j-1
	jp=j+1
	if (ival(i,jm).eq.0) then
		call A_0(w,i,j,ival,idel,a1,a2,a3,a)
		return
	end if
	if (ival(im,j).eq.0) then
		call A0_(w,i,j,ival,idel,a1,a2,a3,a)
		return
	end if
	if (ival(ip,j).eq.0) then
		call AN_(w,i,j,ival,idel,a1,a2,a3,a)
		return
	end if
	if (ival(i,jp).eq.0) then
		call A_N(w,i,j,ival,idel,a1,a2,a3,a)
		return
	end if
	if (ival(im,jm).eq.0) then
		call A_0(0.5d0*w,i,j,ival,idel,a1,a2,a3,a)
		call A0_(0.5d0*w,i,j,ival,idel,a1,a2,a3,a)
		return
	end if
	if (ival(ip,jm).eq.0) then
		call A_0(0.5d0*w,i,j,ival,idel,a1,a2,a3,a)
		call AN_(0.5d0*w,i,j,ival,idel,a1,a2,a3,a)
		return
	end if
	if (ival(im,jp).eq.0) then
		call A_N(0.5d0*w,i,j,ival,idel,a1,a2,a3,a)
		call A0_(0.5d0*w,i,j,ival,idel,a1,a2,a3,a)
		return
	end if
	if (ival(ip,jp).eq.0) then
		call A_N(0.5d0*w,i,j,ival,idel,a1,a2,a3,a)
		call AN_(0.5d0*w,i,j,ival,idel,a1,a2,a3,a)
		return
	end if
	call CORE(w,i,j,ival,idel,im,ip,jm,jp,DX,DY,a1,a2,a3,a)
	return
	end
c
c
c
c
	SUBROUTINE AN_(w,i,j,ival,idel,a1,a2,a3,a)
c
	implicit real*8 (a-h,o-z)
	dimension ival(-1:2,-1:2),idel(-1:2,-1:2),
	1	a1(3,0:1,0:1),a2(3,0:1,0:1),a3(3,0:1,0:1),
	2	a(3,48),DX(0:1,-1:1),DY(0:1,-1:1)
	data DX/0.0,-1.0,	1.0,1.0,	0.0,0.0/
	data DY/0.0,-0.5,	1.0,0.0,	0.0,0.5/
c
	im=i-1
	ip=i
	jm=j-1
	jp=j+1
	if (ival(i,jm).eq.0) then
		call AN0(w,i,j,ival,idel,a1,a2,a3,a)
		return
	end if
	if (ival(i,jp).eq.0) then
		call ANN(w,i,j,ival,idel,a1,a2,a3,a)
		return
	end if
	if (ival(im,jm).eq.0) then
		call AN0(w,i,j,ival,idel,a1,a2,a3,a)
		return
	end if
	if (ival(im,jp).eq.0) then
		call ANN(w,i,j,ival,idel,a1,a2,a3,a)
		return
	end if
	call CORE(w,i,j,ival,idel,im,ip,jm,jp,DX,DY,a1,a2,a3,a)
	return
	end
c
c
c
c
	SUBROUTINE A0N(w,i,j,ival,idel,a1,a2,a3,a)
c
	implicit real*8 (a-h,o-z)
	dimension ival(-1:2,-1:2),idel(-1:2,-1:2),
	1	a1(3,0:1,0:1),a2(3,0:1,0:1),a3(3,0:1,0:1),
	2	a(3,48),DX(0:1,-1:1),DY(0:1,-1:1)
	data DX/0.0,0.0,	1.0,-1.0,	0.0,1.0/
	data DY/0.0,-1.0,	1.0,1.0,	0.0,0.0/
c
	im=i
	ip=i+1
	jm=j-1
	jp=j
	call CORE(w,i,j,ival,idel,im,ip,jm,jp,DX,DY,a1,a2,a3,a)
	return
	end
c
c
c
c
	SUBROUTINE A_N(w,i,j,ival,idel,a1,a2,a3,a)
c
	implicit real*8 (a-h,o-z)
	dimension ival(-1:2,-1:2),idel(-1:2,-1:2),
	1	a1(3,0:1,0:1),a2(3,0:1,0:1),a3(3,0:1,0:1),
	2	a(3,48),DX(0:1,-1:1),DY(0:1,-1:1)
	data DX/0.0,-0.5,	1.0,0.0,	0.0,0.5/
	data DY/0.0,-1.0,	1.0,1.0,	0.0,0.0/
c
	im=i-1
	ip=i+1
	jm=j-1
	jp=j
	if (ival(im,j).eq.0) then
		call A0N(w,i,j,ival,idel,a1,a2,a3,a)
		return
	end if
	if (ival(ip,j).eq.0) then
		call ANN(w,i,j,ival,idel,a1,a2,a3,a)
		return
	end if
	if (ival(im,jm).eq.0) then
		call A0N(w,i,j,ival,idel,a1,a2,a3,a)
		return
	end if
	if (ival(ip,jm).eq.0) then
		call ANN(w,i,j,ival,idel,a1,a2,a3,a)
		return
	end if
	call CORE(w,i,j,ival,idel,im,ip,jm,jp,DX,DY,a1,a2,a3,a)
	return
	end
c
c
c
c
	SUBROUTINE ANN(w,i,j,ival,idel,a1,a2,a3,a)
c
	implicit real*8 (a-h,o-z)
	dimension ival(-1:2,-1:2),idel(-1:2,-1:2),
	1	a1(3,0:1,0:1),a2(3,0:1,0:1),a3(3,0:1,0:1),
	2	a(3,48),DX(0:1,-1:1),DY(0:1,-1:1)
	data DX/0.0,-1.0,	1.0,1.0,	0.0,0.0/
	data DY/0.0,-1.0,	1.0,1.0,	0.0,0.0/
c
	im=i-1
	ip=i
	jm=j-1
	jp=j
	call CORE(w,i,j,ival,idel,im,ip,jm,jp,DX,DY,a1,a2,a3,a)
	return
	end
c
c
c
c
	SUBROUTINE CORE(w,i,j,ival,idel,im,ip,jm,jp,DX,DY,a1,a2,a3,a)
c
	implicit real*8 (a-h,o-z)
	dimension ival(-1:2,-1:2),idel(-1:2,-1:2),
	1	DX(0:1,-1:1),DY(0:1,-1:1),a1(3,0:1,0:1),a2(3,0:1,0:1),
	2	a3(3,0:1,0:1),a(3,48)
c
	jdmax=idel(i,j)/2
	idmax=mod(idel(i,j),2)
c
c	idel:	0 = no derivs
c		1 = x deriv
c		2 = y deriv
c		3 = both derivs
c
	do ja=jm,jp
	do ia=im,ip
	do jd=0,jdmax
	do id=0,idmax
		coeff=w*DX(id,ia-i)*DY(jd,ja-j)
		if (coeff.ne.0.0d0) then
		do kcomp=1,3
			ix=kcomp+3*(ival(ia,ja)-1)
				iy=1
				a(iy,ix)=a(iy,ix)
	1				+a1(kcomp,id,jd)*coeff
				iy=iy+1
				a(iy,ix)=a(iy,ix)
	1				+a2(kcomp,id,jd)*coeff
				iy=iy+1
				a(iy,ix)=a(iy,ix)
	1				+a3(kcomp,id,jd)*coeff
		end do
		end if
	end do
	end do
	end do
	end do
	return
	end
c
c
c
c
c
	SUBROUTINE COEFFS(x,y,a1,a2,a3,xcap,dxxcap,dyxcap)
c
	implicit real*8 (a-h,o-z)
	dimension a1(3,0:1,0:1,0:1,0:1),a2(3,0:1,0:1,0:1,0:1),
	1	a3(3,0:1,0:1,0:1,0:1),xcap(3),dxxcap(3),dyxcap(3),
	2	f_x(0:1,0:1),d_x(0:1,0:1),f_y(0:1,0:1),d_y(0:1,0:1),
	3	p1(3),p2(3)
c
	x1=x
	x2=x**2
	x3=x**3
	f_x(0,0)=1.0d0-3.0d0*x2+2.0d0*x3
	f_x(0,1)=3.0d0*x2-2.0d0*x3
	f_x(1,0)=x1-2.0d0*x2+x3
	f_x(1,1)=-x2+x3
	d_x(0,0)=-6.0d0*x1+6.0d0*x2
	d_x(0,1)=6.0d0*x1-6.0d0*x2
	d_x(1,0)=1.0d0-4.0d0*x1+3.0d0*x2
	d_x(1,1)=-2.0d0*x1+3.0d0*x2

	x1=y
	x2=y**2
	x3=y**3
	f_y(0,0)=1.0d0-3.0d0*x2+2.0d0*x3
	f_y(0,1)=3.0d0*x2-2.0d0*x3
	f_y(1,0)=x1-2.0d0*x2+x3
	f_y(1,1)=-x2+x3
	d_y(0,0)=-6.0d0*x1+6.0d0*x2
	d_y(0,1)=6.0d0*x1-6.0d0*x2
	d_y(1,0)=1.0d0-4.0d0*x1+3.0d0*x2
	d_y(1,1)=-2.0d0*x1+3.0d0*x2

	p2(3)=dsqrt(xcap(1)**2+xcap(2)**2)
	if (p2(3).ne.0.0d0) then
		p1(1)=-xcap(2)/p2(3)
		p1(2)=xcap(1)/p2(3)
	else
		p1(1)=1.0d0
		p1(2)=0.0d0
	end if
	p1(3)=0.0d0
	p2(1)=-xcap(3)*p1(2)
	p2(2)=xcap(3)*p1(1)
	p1dx=0.0d0
	p1dy=0.0d0
	p2dx=0.0d0
	p2dy=0.0d0
	do k=1,3
		p1dx=p1dx+p1(k)*dxxcap(k)
		p1dy=p1dy+p1(k)*dyxcap(k)
		p2dx=p2dx+p2(k)*dxxcap(k)
		p2dy=p2dy+p2(k)*dyxcap(k)
	end do
	area=p1dx*p2dy-p1dy*p2dx
	do jend=0,1
	do iend=0,1
	do jd=0,1
	do id=0,1
		f=f_x(id,iend)*f_y(jd,jend)
		d1=d_x(id,iend)*f_y(jd,jend)*p2dy
	1		-f_x(id,iend)*d_y(jd,jend)*p2dx
		d2=f_x(id,iend)*d_y(jd,jend)*p1dx
	2		-d_x(id,iend)*f_y(jd,jend)*p1dy
		do k=1,3
			if (p2(3).ne.0.0d0) then
				a1(k,id,jd,iend,jend)=f*p2(k)
				a2(k,id,jd,iend,jend)=-f*p1(k)
			else
				a1(k,id,jd,iend,jend)=0.0d0
				a2(k,id,jd,iend,jend)=0.0d0
			end if
			if (area.ne.0.0d0) then
				a3(k,id,jd,iend,jend)=f*xcap(k)
	1				-0.5d0*(d1*p1(k)+d2*p2(k))/area
			else
				a3(k,id,jd,iend,jend)=0.0d0
			end if
		end do
	end do
	end do
	end do
	end do
	return
	end
c
c
c
c
	subroutine INV(imin,n,nsq,lena,asq,y)
c
	implicit real*8 (a-h,o-z)
	dimension lena(n),asq(nsq),y(n)
c
	do i=imin,n
		if (i.gt.imin) then
			im=i-1
			jmin=max0(imin,i-(lena(i)-lena(im))+1)
			do j=jmin,im
				ij=lena(i)-(i-j)
				y(i)=y(i)-asq(ij)*y(j)
			end do
		end if
		ii=lena(i)
		y(i)=y(i)/asq(ii)
	end do
	return
	end
