c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine constraint_sp(jmt,rmt,n_spin_pola,
     >                      vr,r_mesh,pi4,
     >                      evec,evec_r,b_con,b_basis,
     >i_vdif,h_app_para_mag,h_app_perp_mag,iprpts,iprint,istop)
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      implicit   none
c
c     ****************************************************************
c     this routine calculates the constraint necessary for m_x=m_y=0
c     ****************************************************************
c
!     include    'atom_param.h'
c
      character  sname*32
      character  istop*32
c
      integer    i_vdif
      integer    jmt
      integer    n_spin_pola
      integer    j,i
c
      integer iprpts
      integer    iprint
c
      real*4     h_app_para_mag,h_app_perp_mag
      real*4     evec(3)
      real*4     evec_r(3)
      real*4     evec_const(3)
      real*4     b_basis(3,3)
      real*4     b_con_g(3)
      real*4     vr(iprpts,n_spin_pola)
      real*4     vvr(iprpts),br(iprpts)
      real*4     r_mesh(iprpts)
      real*4     b_con(3)
      real*4     rmt
      real*4     b_xcor
      real*4     pi4
      real*4     c
      real*4     b_mag
      real*4     b_con_mag
c
      real*4     zero
      real*4     half
      real*4     one
      real*4     tol
c
      parameter (sname='constraint')
      parameter (zero=0.0d0)
      parameter (half=0.5d0)
      parameter (one=1.0d0)
      parameter (tol=1.0d-6)
c
c     ================================================================
c     Check to see if B_con=0: if so then set evec_r=evec and return..
c     ================================================================
      b_con_mag=zero
      do i=1,3
         b_con_g(i)=b_con(1)*b_basis(i,1)+b_con(2)*b_basis(i,2)+
     >            b_con(3)*b_basis(i,3)
         b_con_mag=b_con_mag+b_con_g(i)*b_con_g(i)
      end do
      b_con_mag=sqrt(b_con_mag)
c
      if(abs(b_con_mag).lt.tol) then
         do i=1,3
            evec_r(i)=evec(i)
         enddo 
         go to 1000
      else if(abs(b_con_mag).gt.one) then
         write(6,*) 'WARNING: b_con_mag is .gt. one !!'
      endif
c
c     ================================================================
c     set up V and B .................................................
c     ----------------------------------------------------------------
      call v_plus_minus_sp(half,vr(1,2),vr(1,1),vvr,br,jmt,iprpts)
c     ----------------------------------------------------------------
c     ================================================================
c     calculate the new direction ....................................
c
      b_mag=0.0d0
      do i=1,3
         evec_r(i)=evec(i)+b_con_g(i)
         b_mag=b_mag+evec_r(i)*evec_r(i)
      end do
      b_mag=sqrt(b_mag)
      do i=1,3
         evec_r(i)=evec_r(i)/b_mag
      end do
c
      write(6,*) 'iprint value :', iprint
      if(iprint.ge.0) then
         write(6,'('' CONSTRAINT: vectors:'')')
         write(6,'(29x,'' Direction:'',18x,''Magnitude :'')')
         write(6,'('' B xcor      :'',3f12.5,5x
     >             )')  (evec(i),i=1,3)
         if(abs(b_con_mag).gt.tol) then
            write(6,'('' B constraint:'',3f12.5,5x,
     >                f10.5)')  (b_con_g(i)/b_con_mag,i=1,3),b_con_mag
         else
            write(6,'('' B constraint:'',3f12.5,5x,
     >                f10.5)')  (b_con_g(i),i=1,3),b_con_mag
         end if
         write(6,'('' B sum       :'',3f12.5,5x,
     >             f10.5)') (evec_r(i),i=1,3),b_mag
      end if
c
c     ================================================================
c     generate new B=sqrt(1+b_con_mag^2)*B_xc.........................
c
      c=sqrt(one+b_con_mag*b_con_mag)
      do j=1,jmt
        br(j)=c*br(j)
      end do
c     ================================================================
c     generate new v_up and v_down....................................
c     ----------------------------------------------------------------
      call v_plus_minus_sp(one,vvr,br,vr(1,2),vr(1,1),jmt,iprpts)
c     ----------------------------------------------------------------
c
 1000 if (istop.eq.sname) then
         call fstop(sname)
      else
         return
      endif
c
      end


