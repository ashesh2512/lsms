c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine getgijmad_sp(mynod,madmat,madmatj,gijmad,ntype,num_atoms,
     >                      map,lmax,kkrsz,
     &               atom_position_1,atom_position_2,atom_position_3,
     >                      rsclu,nrsclu,system_bravais,
     >                      clm,cgnt,
     >                      lofk,mofk,
     >                      ilp1,illp,
     >                      nbortyp,
     >                      pi,pi2,pi4,iprint,istop)
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use MPPModule, only : GlobalCollect, bcastMessage
c
      implicit   none
c
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!      include    'mpif.h'
      include    'param_rsp.h'
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      character  istop*32
      character  sname*20
c
      integer    lmx,kkr,ndlm,ndlj
      parameter  (lmx=1)
      parameter  (kkr=(lmx+1)**2)
      parameter  (ndlj=(2*lmx+1)**2)
      parameter  (ndlm=(2*lmx+1)*(lmx+1))
      integer    maxb
      parameter (maxb=5)
c
      integer     mynod
      integer    lmax,kkrsz
      integer    nrsclu
      integer    lofk(ndlj)
      integer    mofk(ndlj)
      integer    nbortyp(nrsclu)
      integer    itype
      integer    iprint
      integer    ir,is,jr
      integer    i, info
      integer    ntype,num_atoms
      integer    m,n,l
      integer    m1,n1,l1
      integer    map(num_atoms)
      integer    counter
c
      real*4     atom_position_1(num_atoms)
      real*4     atom_position_2(num_atoms)
      real*4     atom_position_3(num_atoms)
      real*4 xoff,yoff,zoff
      real*4     rsclu(iprsclu,3)
      real*4     system_bravais(3,3)
      real*4     clm(ndlm)
      real*4     cgnt((lmax+1)*kkrsz*kkrsz)
      real*4     pi
      real*4     pi2
      real*4     pi4
      real*4     plm(ndlm)
      real*4     sinmp(2*lmx+1)
      real*4     cosmp(2*lmx+1)
      real*4     r2
      real*4     ri(3),rj(3)
      real*4     vects(3)
c
      real*4     zero
      real*4     half
      real*4     one
      real*4     two
      real*4     three
      real*4     tol
      real*4     sqrpi4
c
      real*4 madmat(num_atoms)
      complex*8 madmatj(3,num_atoms,2)
      complex*8 gijmad(kkr,kkr,ntype)
      complex*8 gij(kkr,kkr)
      complex*8 gijtmp(kkr)
      complex*8 hfn(2*lmx+1)
      complex*8 dlm(ndlj)
c
      complex*8 ilp1(0:2*lmx)
      complex*8 illp(kkrsz*kkrsz)
c
      complex*8 czero
      complex*8 cone
c
      parameter (sname='getgijmad')
      parameter (tol=2.d-6)
      parameter (zero=0.0d0)
      parameter (half=0.5d0)
      parameter (one=1.0d0)
      parameter (two=2.0d0)
      parameter (three=3.0d0)
      parameter (czero=(0.0d0,0.0d0))
      parameter (cone=(1.0d0,0.0d0))

      real*4 xp,yp,zp
      real*4 x,y,z,leng2
      leng2(x,y,z)=x*x+y*y+z*z
c
c  ********TESTING**********
      if(iprint.ge.1) then
	do ir=1,num_atoms
	  if(abs(madmatj(2,ir,1)).gt.1.d-10.or.
     &      abs(madmatj(3,ir,1)).gt.1.d-10)
     &      write(6,'(i3,1p4e15.6)') ir,(madmatj(jr,ir,1),jr=2,3)
	enddo
      endif
c  ********TESTING**********

      sqrpi4=sqrt(pi4)
      call zeroout_sp(gijmad,2*kkr*kkr*ntype)
      xoff=atom_position_1(mynod+1)
      yoff=atom_position_2(mynod+1)
      zoff=atom_position_3(mynod+1)

      do ir=2,nrsclu
        itype=nbortyp(ir)
         x=rsclu(ir,1)+xoff
         y=rsclu(ir,2)+yoff
         z=rsclu(ir,3)+zoff
         ri(1)=rsclu(1,1)-rsclu(ir,1)
         ri(2)=rsclu(1,2)-rsclu(ir,2)
         ri(3)=rsclu(1,3)-rsclu(ir,3)

c Find the sublattice this site belongs to
	do is=1,num_atoms
	  if(map(is).eq.itype) then
		 xp=x-atom_position_1(is)
		 yp=y-atom_position_2(is)
		 zp=z-atom_position_3(is)
	       do m1=-maxb,maxb
	       do n1=-maxb,maxb
	       do l1=-maxb,maxb
		 do i=1,3
                 rj(i)=m1*system_bravais(i,1)+n1*system_bravais(i,2)+
     &                 l1*system_bravais(i,3)
		 enddo
		 vects(1)=rj(1)+xp
		 vects(2)=rj(2)+yp
		 vects(3)=rj(3)+zp
	    if(leng2(vects(1),vects(2),vects(3)).le.tol*tol) goto 3
	       enddo ! l1
	       enddo ! n1
	       enddo ! m1
	  endif  ! map
	enddo  ! is
c If we are here something is wrong
	write(6,'(''GETGIJMAD:: error'')')
	stop'getgijmad'
  3         continue
c              =======================================================
c              g(Rij) calculation.....................................
c              -------------------------------------------------------
               call makegij_sp(lmx,kkr,lmx,kkr,
     >                      lmax,kkrsz,ndlj,ndlm,
     >                      czero,ri,sinmp,cosmp,
     >                      clm,plm,cgnt,lofk,mofk,
     >                      ilp1,illp,
     >                      hfn,dlm,gij,
     >                      pi4,iprint,istop)
         madmatj(1,is,1)=madmatj(1,is,1)+gij(1,1)*sqrpi4
c factor of three is (2*l+1)!!
         madmatj(2,is,1)=madmatj(2,is,1)+gij(3,1)
     &       *sqrpi4/three
         madmatj(3,is,1)=madmatj(3,is,1)+gij(4,1)
     &       *sqrpi4/three

      enddo   ! ir

      if(iprint.ge.1) then
	do ir=1,num_atoms
	  if(abs(madmatj(2,ir,1)).gt.1.d-10.or.
     &      abs(madmatj(3,ir,1)).gt.1.d-10)
     &      write(6,'(i3,1p4e15.6)') ir,(madmatj(jr,ir,1),jr=2,3)
	enddo
      endif

c the dipole moments have directions, which may be different for
c different equivalent atoms. We need to find the rotation that rotates
c the atom to one on a real node. We don't really need the rotation.
c Just find the atom that has the same distance to the origin and
c the center of the cluster (within a multiple of the unit vectors).

c  r2 is the distance of the atom ir to the origin
c
c gijmad(*,1,i) is the madelung matrix that generates the dipole
c potential on site i from all other monopoles. gijmad(1,*,i) is
c the madelung matrix that generates the monopole potential on
c site i from all other dipoles.
      do ir=1,num_atoms
	itype=map(ir)
        if (mynod+1 .eq. itype) then
           madmatj(:,:,2) = madmatj(:,:,1)
        end if
!        call mpi_bcast(madmatj(1,1,2), 3*num_atoms, mpi_c_float_complex,
!     &                 itype-1, mpi_comm_world, info)
        call bcastMessage(madmatj(:,1,2),3*num_atoms,  itype-1)
	if(ir.ne.itype) madmat(itype)=madmat(itype)+madmat(ir)
         ri(1)=atom_position_1(ir)-atom_position_1(mynod+1)
         ri(2)=atom_position_2(ir)-atom_position_2(mynod+1)
         ri(3)=atom_position_3(ir)-atom_position_3(mynod+1)
	 r2=leng2(ri(1),ri(2),ri(3))

c factor of two is from Rydberg
	  gijmad(1,1,itype)=gijmad(1,1,itype)+two*madmatj(1,ir,1)
	  gijmad(2,1,itype)=gijmad(2,1,itype)-
     &         two*conjg(madmatj(3,ir,1))
	  gijmad(3,1,itype)=gijmad(3,1,itype)+two*madmatj(2,ir,1)
	  gijmad(4,1,itype)=gijmad(4,1,itype)+two*madmatj(3,ir,1)
	       counter=0
	       call zeroout_sp(gijtmp,2*kkr)
	do jr=1,num_atoms
	  if(map(jr).eq.mynod+1) then
	    xp=atom_position_1(jr)-atom_position_1(itype)
	    yp=atom_position_2(jr)-atom_position_2(itype)
	    zp=atom_position_3(jr)-atom_position_3(itype)
	       do m1=-maxb,maxb
	       do n1=-maxb,maxb
	       do l1=-maxb,maxb
		 do i=1,3
                 rj(i)=m1*system_bravais(i,1)+n1*system_bravais(i,2)+
     &                 l1*system_bravais(i,3)
		 enddo
		 vects(1)=xp+rj(1)
		 vects(2)=yp+rj(2)
		 vects(3)=zp+rj(3)
	    if(abs(r2-leng2(vects(1),vects(2),vects(3))).le.tol) then
	  gijtmp(2)=gijtmp(2)-two*madmatj(3,jr,2)
	  gijtmp(3)=gijtmp(3)+two*conjg(madmatj(2,jr,2))
	  gijtmp(4)=gijtmp(4)+two*conjg(madmatj(3,jr,2))
		 counter=counter+1
	    endif
	       enddo ! l1
	       enddo ! n1
	       enddo ! m1
	  endif ! map
	enddo  ! jr

          if(counter.eq.0) then
	    if(iprint.ge.0) then
	      write(6,'(''getgijmad: no rotation for ir,itype='',3i4)')
     &              ir,itype
	      write(6,'(''getgijmad:: atom_position='',3f10.5)')
     &          atom_position_1(ir),atom_position_2(ir),
     &          atom_position_3(ir),ri
	    endif
	  counter=1
	  gijtmp(2)=gijtmp(2)-two*madmatj(3,ir,1)
	  gijtmp(3)=gijtmp(3)+two*conjg(madmatj(2,ir,1))
	  gijtmp(4)=gijtmp(4)+two*conjg(madmatj(3,ir,1))
	  endif
	 do l=2,kkr
	   gijmad(1,l,itype)=gijmad(1,l,itype)+gijtmp(l)/counter
	 enddo
199      continue
      enddo  !ir

      if(iprint.ge.1) then
	do i=1,ntype
	  write(6,'(''getgijmad:: gijmad: itype='',i4)') i
	  call wrtmtx_sp(gijmad(1,1,i),kkr,istop)
	enddo
      endif
c
c     ================================================================
      if(istop.eq.sname) then
         call fstop(sname)
      else
         return
      endif
c
      end
