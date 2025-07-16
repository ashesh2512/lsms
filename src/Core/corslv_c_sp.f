
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine corslv_sp(jmt,jws,last,last2,h,r,rv,z,
     >                  numc,nc,lc,kc,ecore,
     >                  corden,semden,
     >                  ecorv,esemv,
     >                  n_spin_pola,nrelc,ndeepz,
     >                  iprpts,
     >                  iprint,istop)
c     ================================================================
c
      implicit   none
c
      character  sname*32
      character  istop*32
c
!     include    'atom_param.h'
c
      integer    iprpts
      integer    numc
      integer    nc(numc)
      integer    lc(numc)
      integer    kc(numc)
      integer    jmt
      integer    jws
      integer    n_spin_pola
      integer    nrelc
      integer    ndeepz
      integer    last
      integer    last2
      integer    nitmax
      integer    i
      integer    j
      integer    iter
      integer    iprint
      integer    ndeep
      integer    kappa
      integer    ipdeq

      integer    nnorm
c
      real*4     ecore(numc)
      real*4     r(iprpts),rtmp(0:iprpts)
      real*4     rv(iprpts)
      real*4     corden(iprpts)
      real*4     semden(iprpts)
      real*4     f(0:iprpts)
      real*4     h
      real*4     z
      real*4     zero
      real*4     one
      real*4     two
      real*4     four
      real*4     ten
      real*4     onh
      real*4     tol
      real*4     fnstrc
      real*4     c
      real*4     esemv
      real*4     ecorv
      real*4     fac1
      real*4 gnrm,qmp(0:iprpts)
c
      parameter  (ipdeq=5)
c
      parameter  (zero=0.0d+00)
      parameter  (one=1.0d+00)
      parameter  (two=2.0d+00)
      parameter  (four=4.0d+00)
      parameter  (ten=10.0d+00)
      parameter  (onh=137.0359895d+00)
      parameter  (tol=1.0d-5) ! TOLERANCE FOR SP
      parameter  (sname='corslv')
!     integer solver
!     parameter (solver=1)
      integer nr
      parameter (nr=4)
      real*4 potc(nr,iprpts)
c
c     ****************************************************************
c     this routine drives the radial equations solution for the
c     core states. it must be called for each component(sublattice)
c     separately.
c
c     calls:  deepst(outws,inws)
c             semcor(outws,hank)
c
c     r= radial grid, rv=potential times r
c     ecore= core levels guesses: to be updated
c     g,f upper and lower components times r
c     corden= core charge density
c     nc= principal q. nos., lc= orbital q. nos., kc= kappa q.nos.
c     numc= no. of core states, z= atomic no., jmt=muffin-tin index
c     ndeepz= number of deep core electrons based on zcor
c     jws=top r index
c     ****************************************************************
c
      if (numc.le.0) then
         esemv=zero
         ecorv=zero
         return
      endif
      nitmax=50
!     if(solver.eq.1) then
      fnstrc=one/onh
      c=two/fnstrc
      c=c*ten**nrelc
!     else
!     call fitpot_sp(r,rv,potc,nr,jmt,last)
!     endif
c
c     ================================================================
c     call deepst_sp for each core state.................................
c     ================================================================
      esemv=zero
      ecorv=zero
      ndeep=0
      do i=1,numc
         nnorm = last
         fac1=(3-n_spin_pola)*abs(kc(i))
         ndeep=ndeep+fac1
! if(solver.eq.1) then
c        -------------------------------------------------------------
         call deepst_sp(nc(i),lc(i),kc(i),ecore(i),
     >               rv,r,f(1),h,z,c,nitmax,tol,jws,last,iter,
     >        iprpts,ipdeq)
         if(iter .lt. 0) then
            call fstop('deepst')
         end if
c        -------------------------------------------------------------
c
c        =============================================================
c        when core state energy >-10 then treat as semi-core..........
c        =============================================================
c08/04/95if( ecore(i) .ge. -ten ) then
         if(ndeep.gt.ndeepz)then
            nnorm = last2
c           ----------------------------------------------------------
            call semcst_sp(nc(i),lc(i),kc(i),ecore(i),
     >                    rv,r,f(1),h,z,c,nitmax,tol,jmt,jws,last2,iter,
     >           iprpts,ipdeq)
            if(iter .lt. 0) then
               call fstop('semcst')
            end if
c           ----------------------------------------------------------
         endif
! else
! if(nrelc.eq.0) then
!   kappa=kc(i)
! else
!   kappa=lc(i)
! endif
! call srcore_sp(nrelc.eq.0,ecore(i),f(1),nc(i),kappa,potc,nr,
!    &          r,nitmax,tol,last2,1.d0)
! endif
c     ================================================================
c     normalize the wavefunctions
c     ================================================================
      f(0)=0.d0
      rtmp(0)=0.d0
! meis 4Sep20      do j=1,last2
      do j=1,nnorm
	rtmp(j)=sqrt(r(j))
	f(j)=f(j)/r(j)
      enddo
c     ----------------------------------------------------------------
! meis 4Sep20      call newint_sp(last2+1,rtmp,f,qmp,3)
      call newint_sp(nnorm+1,rtmp,f,qmp,3)
c     ----------------------------------------------------------------
!      gnrm=1.d0/(two*qmp(last2))
      gnrm=1.d0/(two*qmp(nnorm))
      do j=1,last2
         f(j)=f(j)*gnrm*r(j)
      enddo    
         if(ndeep.gt.ndeepz)then
            do j=1,last2
               semden(j)=semden(j) + fac1*f(j)
            enddo
            esemv=esemv+ecore(i)*fac1
	 else
            do j=1,last
! meis 4Sep20             do j=1,last2
               corden(j)= corden(j) + fac1*f(j)
            enddo
            ecorv=ecorv+ecore(i)*fac1
	 endif
      enddo
c
c     ================================================================
c     Major printout: core eigenvalues................................
c     ================================================================
      if(iprint.ge.0) then
         write(6,'(/,'' CORSLV:: Eigenvalues:'',t25,''n'',t30,
     >               ''l'',t35,''k'',t48,''energy'')')
         do i=1,numc
            write(6,'(t23,i3,t28,i3,t33,i3,t41,d20.13)') 
     >      nc(i),lc(i),kc(i),ecore(i)
         enddo
      endif
      if(iprint.ge.0) then
         write(6,'(10x,''Energy     core:'',t40,''='',d16.8)')ecorv
         write(6,'(10x,''Energy semicore:'',t40,''='',d16.8)')esemv
      endif
c     ***************************************************************
      if(iprint.ge.1) then
      write(6,'(/,'' CORSOLV::'')')
         do j=1,last,50
            write(6,'('' j,corden,semden: '',i5,3d20.12)')
     >      j,corden(j),semden(j)
         enddo
      endif
c     ***************************************************************
c
c    ==================================================================
      if(istop.eq.sname) then
      call fstop(sname)
      endif
c
      return
      end