c
c     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine zeroout_sp(x,nx)
c     =================================================================
c
      implicit   none
c
      integer    nx,n
c
      real*4     x(nx)
      real*4     zero
c
      parameter  (zero=0.0)
c
      do n=1,nx
         x(n)=zero
      enddo
c
      return
      end

c
c     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine czeroout_sp(x,nx)
c     =================================================================
c
      implicit   none
c
      integer    nx,n
c
      complex*8     x(nx)
      complex*8     zero
c
      parameter  (zero=(0.0,0.0))
c
      do n=1,nx
         x(n)=zero
      enddo
c
      return
      end

