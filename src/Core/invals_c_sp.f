c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine invals_sp(rg,rf,r,slp,gam,v,z,tol,
     >                  drg,drf,en,c,dk,iprpts,ipdeq)
c     ================================================================
c
c     ****************************************************************
c     initial values for outward integration..........................
c     modified desclaux code.................bg.....june 1990.........
c
c     variables explanation: .........................................
c
c               rg:        big component times r;
c               rf:        small component times r;
c               r:         radial log. grid; 
c               slp:       slope at the origine; 
c               gam:       1st term power;
c               v:         potential at the first point; 
c               z:         atomic number;
c               tol:       check for convergence; 
c               drg,drf:   derivatives times r;
c               c:         speed of light in rydbergs; 
c               en:        energy in rydbergs;
c               dk:        spin angular quantum number;
c     ****************************************************************
      implicit   none
c
!     include    'atom_param.h'
c
      character  sname*20

      integer iprpts,ipdeq
c
      integer    j
      integer    m
      integer    k
c
      real*4     slp
      real*4     gam
      real*4     v
      real*4     z
      real*4     tol
      real*4     en
      real*4     c
      real*4     dk
      real*4     rg(iprpts)
      real*4     rf(iprpts)
      real*4     r(iprpts)
      real*4     drg(2*ipdeq)
      real*4     drf(2*ipdeq)
      real*4     zalfa
      real*4     term1
      real*4     term2
      real*4     pow1
      real*4     pow2
      real*4     dm
      real*4     ratfg
      real*4     sum
      real*4     rfr
      real*4     rgr
      real*4     zero
      real*4     one
      real*4     two
c
      parameter  (zero=0.d0)
      parameter  (one=1.d0)
      parameter  (two=2.d0)
      parameter  (sname='invals_sp')
c

!     write(*,*) "Entering INVALS"
c     ===============================================================
      zalfa=two*z/c
      term1=-zalfa
      term2=zalfa/r(1)+(v-en)/c
c
      do j=1,2*ipdeq
         rg(j)=zero
         rf(j)=zero
      enddo
c
      if( dk .le. zero ) then
         ratfg=(dk-gam)/zalfa
      else
         ratfg=zalfa/(dk+gam)
      endif
c
      rf(2*ipdeq)=slp
      rg(2*ipdeq)=ratfg*slp
c
      do j=1,ipdeq
         rg(j)=rg(2*ipdeq)
         rf(j)=rf(2*ipdeq)
         drg(j)=rg(j)*gam
         drf(j)=rf(j)*gam
      enddo   
      m=1
c
 30   dm=m+gam
      sum=dm*dm-dk*dk+term1*term1
      rfr=(c-term2)*rf(m+2*ipdeq-1)
      rgr=term2*rg(m+2*ipdeq-1)
      pow1=((dm-dk)*rfr-term1*rgr)/sum
      pow2=((dm+dk)*rgr+term1*rfr)/sum
      k=-1
c
      do j=1,ipdeq
         rgr=r(j)**m
         rfr=pow2*rgr
         rgr=pow1*rgr
         if( m.ne.1 .and. abs(rgr/rg(j)).le.tol .and.
     >       abs(rfr/rf(j)).le.tol ) k=1
         rg(j)=rg(j)+rgr
         rf(j)=rf(j)+rfr
         drg(j)=drg(j)+rgr*dm
         drf(j)=drf(j)+rfr*dm
      enddo     
c
      if( k .eq. 1 ) return
c
      rg(m+2*ipdeq)=pow1
      rf(m+2*ipdeq)=pow2
      m=m+1
c
      if( m .le. 20 ) go to 30
c
c     ==============================================================
      write(6,'('' INVALS:: no convergence in small r expansion'')')
      call fstop(sname)
c
      return
      end
