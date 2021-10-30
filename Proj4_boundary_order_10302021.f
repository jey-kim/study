      program boundary_order

      integer xstep,ystep,nplate,xpoints,ypoints,boundary_points,x,y,num
      integer total_points
      double precision rlat,rlon
 
      open (unit=10, file="geometry_kin.dat")
      open (unit=12, file="boundary_order.dat")
      open (unit=14, file="coordinates.dat")

      read(10,*)xstep,ystep,nplate
      
c     xstep = total knotpoints on x axis
c     ystep = total knotpoints on y axis


      xpoints = xstep + 1
      ypoints = ystep + 1 
   
      boundary_points = xpoints*2 + (ypoints-2)*2 
c     boundary_points = total number of points on boundary

      x = xstep
      y = ystep
      num = 0
      write (12,77)boundary_points-1

c fix the y, and change the x (from the top right to top left)
      do i=1,x
      num = num + 1
      xstep = xstep - 1
      ystep = ystep
      write (12,88)num,ystep,xstep
      enddo

c fix the x, and change the y (from the top left to the bottom left)
      do i=1,y
      num = num + 1
      ystep = ystep - 1
      write (12,88)num,ystep,xstep
      enddo

c fix the y, and change the x (from the bottom left to the bottom right)
      do i=1,x
      num = num + 1
      xstep = xstep + 1
      write (12,88)num,ystep,xstep
      enddo

c fix the x and change the y (from the bottom right to the top right)
      do i=1,y
      num = num + 1 
      ystep = ystep + 1
      write (12,88)num,ystep,xstep
      enddo

77    format(1x,I4)
88    format(1x,I4,I4,I4)


      
c extract coordinates from geometry_kin.dat (geometry for kinematic) 

      write(14,*)ypoints,xpoints
      write(14,*)0
      write(14,*)0,3

      total_points = xpoints*ypoints
      
      do i=1,total_points
      read(10,*)
      read(10,*)rlat, rlon
      write(14,99)rlon, rlat
      enddo

99    format(2x, f11.5, 2x,f10.5)

      stop
      end
