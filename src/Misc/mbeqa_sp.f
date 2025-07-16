c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mbeqa_sp(a,b,n)
c     ================================================================
c
      integer    i,n
c
      real*4     a(n)
      real*4     b(n)
c
      do i=1,n
         b(i)=a(i)
      enddo
c
      return
      end