      subroutine gafill_sp(iplmax_ext)
c
!     implicit real*4 (a-h,o-z)
      implicit none
c
!     include 'atom_param.h'
      integer iplmax_ext
      include 'gfill_sp.h'
!     parameter ( lamp=1,lammp=(lamp+1)*(lamp+1) )
c
!     complex*8 gacoeff(lmmaxp,lmmaxp,lammp)
!     complex*8 rgacoeff(kmymaxp,kmymaxp,lammp)
c
!     common/ggaunt/gacoeff
!     common/rggaunt/rgacoeff

      character*32 sname
      character*32 istop
      parameter(sname='gafill')


      integer i,l,m
      integer lam,lp,lmp
      integer lpmax,lpmin
      integer lm,mp
      integer nu
      integer lmmaxp
      integer kkr1

      real*4 gaunt
      real*4 plmg((2*iplmax+1)*(2*iplmax+2)/2,2*iplmax+1)
      real*4 clm((2*iplmax+1)*(2*iplmax+2)/2)
      real*4 tg(2*(2*iplmax+1))
      real*4 wg(2*(2*iplmax+1))
      real*4 sqfpi,fac

c
!     set up things for the gaunt coeff:s
!      call getclm(2*iplmax,clm)
! meis: changed to normalized associated Legendre functions
      call ylm_coefficients_sp(2*iplmax,clm)
      call gauss_legendre_points_sp(-1.d0,1.d0,tg,wg,2*(2*iplmax+1))
      do i=1,2*iplmax+1
! meis: changed to normalized associated Legendre functions
        call plm_normalized_sp(2*iplmax,tg(i),plmg(1,i))
      end do


      lmmaxp =(iplmax+1)*(iplmax+1)

      sqfpi=dsqrt(16.d0*datan(1.d0))
      i=0
      do lam=0,lamp
      fac=sqfpi/(2*lam+1)
      do nu=-lam,lam
        i=i+1
c
        call zeroout_sp(gacoeff(1,1,i),2*lmmaxp*lmmaxp)
        call zeroout_sp(rgacoeff(1,1,i),2*kmymaxp*kmymaxp)
c
        do l=0,iplmax
        do m=-l,l
          lm=l*l+l+m+1
          mp=m-nu
          lpmin=iabs(l-lam)
          lpmax=l+lam
          do 1 lp=lpmin,lpmax,2
            if(lp.gt.iplmax) goto 1
            if(iabs(mp).gt.lp) goto 1
            lmp=lp*lp+lp+mp+1
!           gacoeff(lm,lmp,i)=fac*skkr_gaunt(l,lam,lp,m,nu,mp)
            gacoeff(lm,lmp,i)=fac*gaunt(lam,nu,l,m,lp,mp,
     >                         2*iplmax,2*iplmax+1,wg,plmg,clm)
    1     continue
        end do
        end do
c       call outmat_sp(gacoeff(1,1,i),lmmaxp,lmmaxp,lmmaxp,6)
!       call relmtrx_sp(gacoeff(1,1,i),rgacoeff(1,1,i),iplmax)
        kkr1=(iplmax+1)*(iplmax+1)
        call relmtrx_sp(gacoeff(1,1,i),rgacoeff(1,1,i),kkr1,kkr1)
c       call outmat_sp(rgacoeff(1,1,i),kmymaxp,kmymaxp,kmymaxp,6)
c
      end do
      end do

c
      return
      end