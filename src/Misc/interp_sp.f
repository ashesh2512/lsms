c
c
c
      subroutine interp_sp(r,f,nr,rs,ps,dps,deriv)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c----- this routine interpolates functn f on grid r
c      returns value ps and its derivative dps 
c
      implicit real*4 (a-h,o-z)
      logical deriv
      real*4 r(nr),f(nr),c(4)
c
      ip1=2
      do i=ip1,nr-1
      if (rs.gt.r(i)) ip1=i+1
      enddo
      i=ip1-1
      call fit_sp(r,f,nr,i,c)
c
      h1=rs-r(i)
      h2=r(ip1)-rs
      c2=1.d0+c(4)*h1*h2
      ps=((c(3)-c(2))*h1+c(1))/c2
      dps=c(2)+(c(3)-c(2)+ps*c(4)*(h1-h2))/c2
      ps=ps+f(i)-c(1)+c(2)*h1
      
      return
      end
