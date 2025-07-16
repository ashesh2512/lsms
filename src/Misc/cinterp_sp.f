c
      subroutine cinterp_sp(r,f,nr,rs,ps,dps,work)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c----- this routine interpolates functn f on grid r
c      returns value ps and its derivative dps 
c
      implicit real*4 (a-h,o-z)
      logical deriv
      real*4 r(nr),c(4),work(nr),psr,psi,dpsr,dpsi
      complex*8 f(nr),ps,dps
c
      ip1=2
      do i=ip1,nr-1
      if (rs.gt.r(i)) ip1=i+1
      enddo

c Real part
      do i=1,nr
	work(i)=real(f(i))
      enddo
      i=ip1-1
      call fit_sp(r,work,nr,i,c)
c
      h1=rs-r(i)
      h2=r(ip1)-rs
      c2=1.d0+c(4)*h1*h2
      psr=((c(3)-c(2))*h1+c(1))/c2
      dpsr=c(2)+(c(3)-c(2)+psr*c(4)*(h1-h2))/c2
      psr=psr+real(f(i))-c(1)+c(2)*h1

c Imaginary part
      do i=1,nr
	work(i)=imag(f(i))
      enddo
      i=ip1-1
      call fit_sp(r,work,nr,i,c)
c
      h1=rs-r(i)
      h2=r(ip1)-rs
      c2=1.d0+c(4)*h1*h2
      psi=((c(3)-c(2))*h1+c(1))/c2
      dpsi=c(2)+(c(3)-c(2)+psi*c(4)*(h1-h2))/c2
      psi=psi+imag(f(i))-c(1)+c(2)*h1

      ps=dcmplx(psr,psi)
      dps=dcmplx(dpsr,dpsi)
      
      return
      end