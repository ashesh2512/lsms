! single_site_tmat provides an interface for the single site
! scattering routines. 

      subroutine single_site_tmat_sp(nrel_rel,n_spin_cant,is,
     >                       n_spin_pola,
     >                       mtasa,rws,
     >                       nrelv,clight,lmax,kkrsz,
     >                       energy,prel,pnrel,
     >                       vr,h,jmt,jws,r_mesh,
     >                       tmat_l,tmat_g,matom,
     >                       zlr,jlr,
     >                       gz,fz,gj,fj,nuz,indz,
     >                       ubr,ubrd,dmat,dmatp,
     >                       r_sph,iprint,istop)

! tmat_l returns the spin-diagonal blocks of the tmatrix in the local frame
!        in the non-relativistic case only. Its value in the
!        relativistic case is undefined.
! tmat_g returns the tmatrix in the global frame of reference.

      implicit none

      include    'atom_param.h'

      character  sname*32
      character  istop*32
      character  idpot*10

      integer    nuzp
      parameter  (nuzp=2)

      integer    nrel_rel     ! nrel_rel .ne. 0 : relativistic calc.
      integer    n_spin_cant
      integer    n_spin_pola
      integer    is
      integer    nrelv
      integer    mtasa
      integer    lmax
      integer    kkrsz
      integer    jmt
      integer    jws
      integer    iprint
      integer    nuz(2*kkrsz)
      integer    indz(nuzp,2*kkrsz)
      integer    kmymax
      integer    j
      integer    info
      integer    ipvt(2*kkrsz)
      integer    isq
      integer    isp

      real*4     vr(iprpts,n_spin_pola),v0
      real*4     clight
      real*4     h
      real*4     r_sph
      real*4     rws
      real*4     r_mesh(iprpts)
      real*4     vrr(iprpts)
      real*4     brr(iprpts)
      real*4     boprr(iprpts,2)
      real*4     soscal

      complex*8 matom(lmax+1,n_spin_pola)
      complex*8 tmat_l(kkrsz*kkrsz,n_spin_cant)
      complex*8 tmat_g(kkrsz*n_spin_cant,kkrsz*n_spin_cant)
      complex*8 energy
      complex*8 prel
      complex*8 psq
      complex*8 pnrel
      complex*8 zlr(iprpts*(iplmax+1),ip_spin_cant)
      complex*8 jlr(iprpts*(iplmax+1),ip_spin_cant)
      complex*8 gz(iprpts,nuzp,2*kkrsz)
      complex*8 fz(iprpts,nuzp,2*kkrsz)
      complex*8 gj(iprpts,nuzp,2*kkrsz)
      complex*8 fj(iprpts,nuzp,2*kkrsz)
      complex*8 wbig(4*kkrsz*kkrsz)
      complex*8 ubr(4)
      complex*8 ubrd(4)
      complex*8 dmat(2*kkrsz,2*kkrsz)
      complex*8 dmatp(2*kkrsz,2*kkrsz)
      complex*8 detl

      parameter (sname='single_site_tmat')
! soscal: spin-orbit scaling:
! soscal=1 : Dirac eq with full spin-orbit coupling
!        0 : S.O. completly scaled out
      parameter (soscal=1.d0)
      idpot='x123456789'
c     if(vr(1,1).gt.-1.d-3) idpot='Vacuum    '
c     if(vr(1,1).gt.-1.d-3) v0=vr(1,1)/r_mesh(1)
      v0=0.d0
c     if(vr(1,1).gt.-1.d-3) write(*,*)idpot
!     write(6,*) 'entering ',sname
      kmymax = kkrsz*n_spin_cant
      if(nrel_rel.eq.0) then
!       -------------------------------------------------
!        non relativistic
!       -------------------------------------------------
        do isp=1,n_spin_cant
          isq=max(is,isp)
c         -------------------------------------------------------------
          if(iprint.ge.1) then
            write(6,*) 'before SINGLE_SCATTERER_NONREL'
            call flush(6)
          end if
          call single_scatterer_nonrel_sp(nrelv,clight,lmax,kkrsz,
     >                     energy,prel,pnrel,
     >                     vr(1,isq),r_mesh,h,jmt,jws,
     >                     tmat_l(1,isp),matom(1,isp),
     >                      zlr(1,isp),jlr(1,isp),
     >                     r_sph,iprint,istop)
c         -------------------------------------------------------------
          if(iprint.ge.1) then
c           ===========================================================
c           write out Single-site t-matrix (local frame)...............
            write(6,'(/,''GETTAU:: t-matrixfor isp='',i5)') isp
c           -----------------------------------------------------------
            call wrtmtx_sp(tmat_l(1,isp),kkrsz,istop)
c           -----------------------------------------------------------
          endif
        enddo

        if(n_spin_cant.eq.2) then
          call trltog_sp(kkrsz,kkrsz,ubr,ubrd,
     >               tmat_l(1,1),tmat_l(1,2),tmat_g)
        else
          call ccopy(kkrsz*kkrsz,tmat_l(1,1),1,tmat_g,1)
        end if

      else
!       -------------------------------------------------
!       relativistic
!       -------------------------------------------------
!       first build the Clebsch - Gordan coefs (should be done better
!       and not here) !!!!
!       write(6,*) sname,1
!       call clebsch_sp  !moved to beginning of zplanint
!       currently the relativistic part works only for ASA:
        if(mtasa.ne.1) then
          write(6,*) '  relativistic calcultions work for ASA only!'
          call fstop(sname)
        end if
!       first we transform the potential from the up-down form
!        to el-mag form
!       write(6,*) sname,2
        do j=1,jws
!         do we need a factor of 1/2 here or not?
          vrr(j) = 0.5d0*(vr(j,1)+vr(j,2))
          if(vr(1,1).gt.-1.d-4)vrr(j)=vrr(j)-1.d-4
          brr(j) = 0.5d0*(vr(j,1)-vr(j,2))
!         vrr(j) =       (vr(j,1)+vr(j,2))
!         brr(j) =       (vr(j,1)-vr(j,2))
          boprr(j,1) = 0.d0
          boprr(j,2) = 0.d0
        end do
!         write(6,*) h,jws,rws
!         write(6,200) (vrr(j),j=1,jws)
!         write(6,200) (brr(j),j=1,jws)
  200   format(4d20.12)
        psq = energy  +energy*energy/(clight*clight)
!       h = 1.200000000000000d-002
!       rws = 2.60957430788315d0
        call single_scatterer_rel_sp(energy,psq,lmax,kmymax,
     >            'x123456789',v0,
     >            vrr,brr,boprr,h,jws,rws,
     >            tmat_g,gz,fz,gj,fj,nuz,indz,1,soscal,
     >            iprint,istop)

!       now we have tinv but we need t:
!       write(6,*) sname,4
!       call cgetrf(kmymax,kmymax,tmat_g,kkrsz,ipvt,info)
!       call cgetri(kmymax,tmat_g,kmymax,ipvt,wbig,kkrsz*kkrsz,info)
        call gjinv_sp(tmat_g,kmymax,kmymax,detl)
!       t now is in the local frame of refference, now we have to
!       transform it into the global frame:
!       write(6,*) sname,5
        call tripmt_sp(dmat,tmat_g,dmatp,kmymax,kmymax,kmymax)

      end if

      if(istop.eq.sname) then
        write(6,*) 'tmat_g:'
        call outmat1_sp(tmat_g,kmymax,kmymax,kmymax,1.d-15,6)
        call gjinv_sp(tmat_g,kmymax,kmymax,detl)
        call replms_sp(wbig,tmat_g,lmax,kmymax)
        write(6,*) 'diag(tmat_lms):'
        do j=1,2*(lmax+1)**2
          write(6,*) j,j,wbig(j)
        end do
        write(6,*) 't^(-1) :'
        call outmat1_sp(tmat_g,kmymax,kmymax,kmymax,1.d-15,6)
        call fstop(sname)
      end if

      end
