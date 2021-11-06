      program boundary_spline_fit

      implicit real*8 (a-h,o-z)
 
      open (unit=10, file="spline_fit_null.dat")
      open (unit=12, file="boundary_spline_fit.dat")

      read(10,*)nx,ny,ntotal
      read(10,*)l
      read(10,*)ll,dummy1,dummy2,dummy3
      read(10,*)ngrid

      write(12,77)nx,ny,ntotal
      write(12,77)117
      write(12,77)ngrid

      do i=1,ngrid
         read(10,*)
         read(10,*)ii,jj,kk
         read(10,*)exx,eyy,exy
         read(10,*)sxx,syy,sxy,cc1,cc2,cc3
         write(12,77)ii,jj,kk
         write(12,88)exx,eyy,exy
         write(12,88)sxx,syy,sxy,cc1,cc2,cc3
      enddo

77    format(1x,I5,1x,I5,1x,I5)
88    format(6e15.5)
      end
