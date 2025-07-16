c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine cal_madelung_matrix_sp(mynod,num_atoms,
     >                               bravais_lattice_in,
     >                               atom_posi_x_in,
     >                               atom_posi_y_in,
     >                               atom_posi_z_in,
     >                               madmat,
     >                               iprint,istop)
c     ================================================================
c
      implicit   none
c
!     include 'atom_param.h'
      include 'madelung.h'
c
      character  istop*32
      character  sname*32
      parameter (sname='cal_madelung_matrix')
c
      integer    mynod
      integer    num_atoms
c
      integer    iprint
      integer    i
      integer    imin
      integer    nrslat
      integer    nknlat
c
      integer    j,k
c
      real*4     bravais_lattice_in(9)
      real*4     atom_posi_x_in(num_atoms)
      real*4     atom_posi_y_in(num_atoms)
      real*4     atom_posi_z_in(num_atoms)
      real*4     bravais_lattice(9)
      real*4     atom_posi_x(num_atoms)
      real*4     atom_posi_y(num_atoms)
      real*4     atom_posi_z(num_atoms)
      real*4     madmat(num_atoms)
c
      real*4     a0
      real*4     dmin
      real*4     dmin2
      real*4     dmax
      real*4     dmax2
      real*4     dx
      real*4     dy
      real*4     dz
      real*4     dis2
      real*4     dis
      real*4     a1
      real*4     a2
      real*4     a3
      real*4     omegbra
      real*4     rslat_x(iprslat_mad)
      real*4     rslat_y(iprslat_mad)
      real*4     rslat_z(iprslat_mad)
      real*4     rslatsq(iprslat_mad)
      real*4     knlat_x(ipknlat_mad)
      real*4     knlat_y(ipknlat_mad)
      real*4     knlat_z(ipknlat_mad)
      real*4     knlatsq(ipknlat_mad)
      real*4     etainv
      real*4     eta
      real*4     etainv0
      real*4     sum
      real*4     one
      real*4     zero
      parameter (one=1.0d0)
      parameter (zero=0.0d0)
c
      parameter (etainv0=0.5d0)
c
!     if(max_atoms.lt.num_atoms) then
!        write(6,'('' CAL_MADELUNG_MATRIX:: max_atoms < num_atoms'',
!    >             2i5)')max_atoms,num_atoms
!        call fstop('cal_madelung_matrix')
!     endif
c     ----------------------------------------------------------------
      call mbeqa_sp(bravais_lattice_in,bravais_lattice,9)
      call mbeqa_sp(atom_posi_x_in,atom_posi_x,num_atoms)
      call mbeqa_sp(atom_posi_y_in,atom_posi_y,num_atoms)
      call mbeqa_sp(atom_posi_z_in,atom_posi_z,num_atoms)
c     ----------------------------------------------------------------
c
      a1=sqrt(bravais_lattice(1)*bravais_lattice(1)+
     >        bravais_lattice(2)*bravais_lattice(2)+
     >        bravais_lattice(3)*bravais_lattice(3) )
      a2=sqrt(bravais_lattice(4)*bravais_lattice(4)+
     >        bravais_lattice(5)*bravais_lattice(5)+
     >        bravais_lattice(6)*bravais_lattice(6) )
      a3=sqrt(bravais_lattice(7)*bravais_lattice(7)+
     >        bravais_lattice(8)*bravais_lattice(8)+
     >        bravais_lattice(9)*bravais_lattice(9) )
      a0=min(a1,a2,a3)
ctest
      etainv=etainv0+0.01*max(a1,a2,a3)/a0
      etainv=etainv*a0
ctest
c     ================================================================
c     change units so that both bravais_lattice and atom_posi_* are in
c     in the units of a0 = 1
c     ----------------------------------------------------------------
c     call dscal_sp(9,one/a0,bravais_lattice,1)
c     call dscal_sp(num_atoms,one/a0,atom_posi_x,1)
c     call dscal_sp(num_atoms,one/a0,atom_posi_y,1)
c     call dscal_sp(num_atoms,one/a0,atom_posi_z,1)
c     ----------------------------------------------------------------
c
c     ================================================================
c     obtain the lattice vectors for the big cell.....................
c     rslat, rslatsq, knlat, and knlatsq are in the units of a0 = 1...
c     ----------------------------------------------------------------
      call getstruc_sp(mynod,bravais_lattice,omegbra,etainv,
     >              iprslat_mad,ipknlat_mad,
     >              rslat_x,rslat_y,rslat_z,rslatsq,nrslat,
     >              knlat_x,knlat_y,knlat_z,knlatsq,nknlat,
     >              iprint,istop)
c     ----------------------------------------------------------------
c
c     ================================================================
c     set up the Madelung matrix......................................
c     ----------------------------------------------------------------
      call madewd_sp(nrslat,nknlat,mynod+1,num_atoms,one,
     >            rslat_x,rslat_y,rslat_z,rslatsq,
     >            knlat_x,knlat_y,knlat_z,knlatsq,
     >            atom_posi_x,atom_posi_y,atom_posi_z,
     >            etainv,omegbra,madmat,
     >            iprint,istop)
c     ----------------------------------------------------------------
c     ================================================================
      if(sname.eq.istop) then
         call fstop(sname)
      else
         return
      endif
      end