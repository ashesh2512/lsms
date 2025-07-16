c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine inwhnk_sp(invp,nmax,rg,rf,rv,r,en,c,drg,drf,
     >                  dk,dm,nws,imm,iprpts,ipdeq)
c     ================================================================
      implicit   none
c
!     include   'atom_param.h'
c
c     ****************************************************************
c     drg,drf:   derivative of dp and dq times r;  
c     c      :   speed of light.
c
c     free solutions ....... riccati-hankel functions 
c     ****************************************************************
c
      integer    iprpts,ipdeq
      integer    invp
      integer    nmax
      integer    nws
      integer    imm
      real*4    kappa
      integer    l
      integer    ll
      integer    lp
      integer    n
c
      real*4     rg(iprpts)
      real*4     rf(iprpts)
      real*4     rv(iprpts)
      real*4     r(iprpts)
      real*4     en
      real*4     c
      real*4     drf(2*ipdeq)
      real*4     drg(2*ipdeq)
      real*4     dk
      real*4     dm
c
      real*4     mass
      real*4     zero
      real*4     cinv
      real*4     c2inv
      real*4     ak
      real*4     sgnk
      real*4     factor
c
      complex*8 bh(10)
      complex*8 dh(10)
      complex*8 x
      complex*8 p
      complex*8 mil
      complex*8 sqrtm1
      complex*8 mi
c
      parameter (sqrtm1=(0.d0,1.d0))
      parameter (mi=(0.d0,-1.d0))
      parameter (zero=0.d0)
c
      cinv=1.d0/c
      c2inv=cinv*cinv 
      mass=(1.d0 + en*c2inv)
      p=sqrt(dcmplx(en*mass,zero) )
c
      kappa=dk
      ak=abs(dk)
      sgnk=dk/ak
      factor= sgnk*sqrt( -en/mass )*cinv
c
c     ================================================================
c     officially this is L+1
c     ================================================================
      l = ak + (1+sgnk)/2
      lp=l - sgnk
      mil=sqrtm1**(l-1)
c
c     ================================================================
c     maximum l needed for this kappa
c     ================================================================
      ll=max(l,lp)
c
      do n=iprpts,invp,-1
         x=p*r(n)
c        -------------------------------------------------------------
         call richnk_sp(ll,x,bh,dh)
c        -------------------------------------------------------------
c
c        =============================================================
c        NOTE: working on log-mesh so need extra factor of r for rdg 
c              and drf.
c        =============================================================
         rg(n)=-mi*mil*bh(l)
         rf(n)=-factor*mil*bh(lp)
      enddo
c
      drg(ipdeq)=-x*(mi*mil*dh(l))
      drf(ipdeq)=-x*factor*mil*dh(lp)
c
      return
      end
