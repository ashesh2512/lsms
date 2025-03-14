/* -*- c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */
#include <cmath>
#include <vector>

#include <hdf5.h>

#include "Complex.hpp"
#include "PhysicalConstants.hpp"
#include "Real.hpp"

// #include "lapack.h"

// #include "Communication/LSMSCommunication.hpp"
#include "Main/LSMSMode.hpp"
#include "Misc/Coeficients.hpp"
#include "Misc/Indices.hpp"
#include "MultipleScattering.hpp"
#include "SingleSite/SingleSiteScattering.hpp"
#include "buildKKRMatrix.hpp"
#include "linearSolvers.hpp"
#include "tau00Postprocess.cpp"
#ifdef USE_NVTX
#include <nvToolsExt.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#else
#ifndef LSMS_DUMMY_OPENMP
#define LSMS_DUMMY_OPENMP
inline int omp_get_max_threads() { return 1; }
inline int omp_get_num_threads() { return 1; }
inline int omp_get_thread_num() { return 0; }
#endif
#endif

#if defined(ACCELERATOR_LIBSCI) || defined(ACCELERATOR_CUDA_C) || \
    defined(ACCELERATOR_HIP)
// #include <cuda_runtime.h>
#include "Accelerator/DeviceStorage.hpp"
// extern void * deviceStorage;
extern DeviceStorage *deviceStorage;
#endif
#include <assert.h>

extern "C" {
void write_kkrmat_(Complex *a, int *n, int *lda, Complex *e);
};

// void buildKKRMatrix_gpu(LSMSSystemParameters &lsms, LocalTypeInfo
// &local,AtomData &atom, Complex energy, Complex prel, int iie, Matrix<Complex>
// &m);
void block_inverse(Matrix<Complex> &a, int *blk_sz, int nblk,
                   Matrix<Complex> &delta, int *ipvt, int *idcol);

// #define SYNTHETIC_MATRIX
// #define WRITE_GIJ

void buildKKRMatrix(LSMSSystemParameters &lsms, LocalTypeInfo &local,
                    AtomData &atom, int ispin, Complex energy, Complex prel,
                    int iie, Matrix<Complex> &m) {
  Real rij[3];

  int lmax = lsms.maxlmax;
  int kkrsz = (lmax + 1) * (lmax + 1);
  int kkrsz_ns = kkrsz * lsms.n_spin_cant;
  int nrmat_ns = lsms.n_spin_cant * atom.nrmat;
  int nrst, ncst;

  Complex *gij = new Complex[kkrsz * kkrsz];
  Real *sinmp = new Real[2 * lmax + 1];
  Real *cosmp = new Real[2 * lmax + 1];
  Real *plm = new Real[AngularMomentumIndices::ndlm];
  Complex *hfn = new Complex[2 * lmax + 1];
  Complex *dlm = new Complex[AngularMomentumIndices::ndlj];
  Complex *bgij = new Complex[4 * kkrsz * kkrsz];
  Complex *tmat_n = new Complex[atom.kkrsz * atom.kkrsz * 4];

  const Complex cmone = -1.0;
  const Complex czero = 0.0;

  Real pi4 = 4.0 * 2.0 * std::asin(1.0);

  for (int i = 0; i < nrmat_ns * nrmat_ns; i++) m[i] = 0.0;
  for (int i = 0; i < nrmat_ns; i++) m(i, i) = 1.0;

  // print first t matrix
  if (lsms.lsmsMode == LSMSMode::liz0 && lsms.rank == 0) {
    printf("first t Matrix:\n");
    for (int i = 0; i < kkrsz_ns; i++)
      for (int j = 0; j < kkrsz_ns; j++) {
        printf("%d %d %.8le %.8le\n", i, j,
               std::real(local.tmatStore(i + j * kkrsz_ns, 0)),
               std::imag(local.tmatStore(i + j * kkrsz_ns, 0)));
      }
    // write the tmatStore
    FILE *of = fopen("tmat.dat", "w");
    fprintf(of, "# tmatStore file for buildkkrmat test:\n");
    fprintf(of, "# line 4: num_store kkrsz Re(energy) Im(energy)\n");
    fprintf(of,
            "# following numstore*kkrsz*kkrsz lines: storeidx i j Re(t_ij) "
            "Im(t_ij)\n");
    fprintf(of, "%4zu %4d %lg %lg\n", local.tmatStoreGlobalIdx.size(), kkrsz_ns,
            std::real(energy), std::imag(energy));
    for (int idx = 0; idx < local.tmatStoreGlobalIdx.size(); idx++) {
      for (int i = 0; i < kkrsz_ns; i++)
        for (int j = 0; j < kkrsz_ns; j++) {
          fprintf(of, "%4d %4d %4d %.8le %.8le\n", idx, i, j,
                  std::real(local.tmatStore(i + j * kkrsz_ns, idx)),
                  std::imag(local.tmatStore(i + j * kkrsz_ns, idx)));
        }
    }

    fclose(of);
    // write LIZ positions for first atom
    //    of=fopen("LIZ_pos.h","w");
    //    fprintf(of,"atom.numLIZ=%d;\n",atom.numLIZ);
    //    fprintf(of,"atom.LIZPos.resize(3,atom.numLIZ);\n");
    //    for(int i=0; i<atom.numLIZ; i++)
    //    {
    //      fprintf(of,"atom.LIZPos(0,%d)=%lg;\n",i,atom.LIZPos(0,i));
    //      fprintf(of,"atom.LIZPos(1,%d)=%lg;\n",i,atom.LIZPos(1,i));
    //      fprintf(of,"atom.LIZPos(2,%d)=%lg;\n",i,atom.LIZPos(2,i));
    //    }
    //    fclose(of);

    of = fopen("LIZ_pos.dat", "w");
    fprintf(of, "# LIZ file for buildkkrmat test:\n");
    fprintf(of, "# line 4: LIZsize\n");
    fprintf(of, "# following LIZsize line: i x_i y_i z_i idx_i lmax_i\n");
    fprintf(of, "%d\n", atom.numLIZ);
    for (int i = 0; i < atom.numLIZ; i++) {
      fprintf(of, "%6d ", i);
      fprintf(of, "%lg ", atom.LIZPos(0, i));
      fprintf(of, "%lg ", atom.LIZPos(1, i));
      fprintf(of, "%lg ", atom.LIZPos(2, i));
      fprintf(of, "%4d ", atom.LIZStoreIdx[i]);
      fprintf(of, "%2d\n", atom.LIZlmax[i]);
    }
    fclose(of);

    exit(0);
  }

#ifdef WRITE_GIJ
  Matrix<Complex> Gij_full(nrmat_ns, nrmat_ns);
  for (int i = 0; i < nrmat_ns * nrmat_ns; i++) Gij_full[i] = 0.0;
#endif

  nrst = 0;
  for (int ir1 = 0; ir1 < atom.numLIZ; ir1++) {
    int kkr1 = (atom.LIZlmax[ir1] + 1) * (atom.LIZlmax[ir1] + 1);
    int kkr1_ns = kkr1 * lsms.n_spin_cant;

    ncst = 0;
    // build local t_mat:
    int im = 0;
    if (lsms.n_spin_pola == lsms.n_spin_cant)  // non polarized or spin canted
    {
      for (int js = 0; js < lsms.n_spin_cant; js++) {
        int jsm = kkrsz * kkrsz_ns * js;
        for (int j = 0; j < kkr1; j++) {
          for (int is = 0; is < lsms.n_spin_cant; is++) {
            int jm = jsm + kkrsz_ns * j + kkrsz * is;
            int one = 1;
            BLAS::zcopy_(&kkr1,
                         &local.tmatStore(iie * local.blkSizeTmatStore + jm,
                                          atom.LIZStoreIdx[ir1]),
                         &one, &tmat_n[im], &one);
            im += kkr1;
          }
        }
      }
    } else {  // spin polarized colinear version for ispin
      int jsm = kkrsz * kkrsz * ispin;  // copy spin up or down?
      for (int j = 0; j < kkr1; j++) {
        int jm = jsm + kkrsz_ns * j;
        int one = 1;
        BLAS::zcopy_(&kkr1,
                     &local.tmatStore(iie * local.blkSizeTmatStore + jm,
                                      atom.LIZStoreIdx[ir1]),
                     &one, &tmat_n[im], &one);
        im += kkr1;
      }
    }
    /*
        if(ir1==1)
        {
          printf("atom.LIZlmax[ir1]=%d\n",atom.LIZlmax[ir1]);
          for(int i=0; i<atom.kkrsz*atom.kkrsz*4; i++)
          {
            Complex t=local.tmatStore(i,atom.LIZStoreIdx[ir1]);
            if(std::norm(t)!=0.0)
            {
              printf("tmat_n[%d]=(%lf,
       %lf)\n",i,std::real(tmat_n[i]),std::imag(tmat_n[i]));
              printf("tmatStore[%d,ir]=(%lf,
       %lf)\n",i,std::real(t),std::imag(t));
            }
          }
        }
    */
    //
    for (int ir2 = 0; ir2 < atom.numLIZ; ir2++) {
      int kkr2 = (atom.LIZlmax[ir2] + 1) * (atom.LIZlmax[ir2] + 1);
      int kkr2_ns = kkr2 * lsms.n_spin_cant;
      if (ir1 != ir2) {
        rij[0] = atom.LIZPos(0, ir1) - atom.LIZPos(0, ir2);
        rij[1] = atom.LIZPos(1, ir1) - atom.LIZPos(1, ir2);
        rij[2] = atom.LIZPos(2, ir1) - atom.LIZPos(2, ir2);

        makegij_(&atom.LIZlmax[ir1], &kkr1, &atom.LIZlmax[ir2], &kkr2,
                 &lsms.maxlmax, &kkrsz, &AngularMomentumIndices::ndlj,
                 &AngularMomentumIndices::ndlm, &prel, &rij[0], sinmp, cosmp,
                 &SphericalHarmonicsCoeficients::clm[0], plm,
                 &GauntCoeficients::cgnt(0, 0, 0), &GauntCoeficients::lmax,
                 &AngularMomentumIndices::lofk[0],
                 &AngularMomentumIndices::mofk[0], &IFactors::ilp1[0],
                 &IFactors::illp(0, 0), hfn, dlm, gij, &pi4,
                 &lsms.global.iprint, lsms.global.istop, 32);
        Complex psq = prel * prel;
        int nrel_rel = 0;
        if (lsms.relativity == full) nrel_rel = 1;
        setgij_(gij, bgij, &kkr1, &kkr1_ns, &kkr2, &kkr2_ns, &lsms.n_spin_cant,
                &nrel_rel, &psq, &energy);

//        if((ir1==1 && ir2==0) || (ir1==10 && ir2==0))
//        {
//          printf("ORIG: ir1=%d, ir2=%d: bgij[0] = %g + %gi; gij[0] = %g +
//          %gi\n",
//                 ir1, ir2, bgij[0].real(), bgij[0].imag(), gij[0].real(),
//                 gij[0].imag());
//          printf("    rij = %g %g %g;  prel=%g + %gi\n", rij[0],  rij[1],
//          rij[2], prel.real(), prel.imag()); printf("    kkr1 = %d; kkr2 = %d;
//          kkrsz = %d\n", kkr1, kkr2, kkrsz);
//        }
#ifdef WRITE_GIJ
        for (int ii = nrst; ii < nrst + kkr1_ns; ii++)
          for (int jj = ncst; jj < ncst + kkr2_ns; jj++)
            Gij_full(ii, jj) = bgij[ii - nrst + (jj - ncst) * kkr2_ns];
#endif
        BLAS::zgemm_("n", "n", &kkr1_ns, &kkr2_ns, &kkr1_ns, &cmone, tmat_n,
                     &kkr1_ns, bgij, &kkr1_ns, &czero, &m(nrst, ncst),
                     &nrmat_ns);
      }
      ncst += kkr2_ns;
    }
    nrst += kkr1_ns;
  }

#ifdef SYNTHETIC_MATRIX
  std::cout << "WARNING: USING SYNTHETIC MATRIX!!!!\n";
  for (int i = 0; i < nrmat_ns; i++)
    for (int j = 0; j < nrmat_ns; j++) {
      m(i, j) = Complex((11.0 * (i + 1) - 3.0 * (j + 1)) / Real(i + j + 2),
                        (5.0 * (i + 1) - 2.0 * (j + 1)) / Real(i + j + 2));
    }
#endif

  /*
  #ifdef WRITE_GIJ
    write_kkrmat_(&Gij_full(0,0),&nrmat_ns,&nrmat_ns, &energy);
  #else
    write_kkrmat_(&m(0,0),&nrmat_ns,&nrmat_ns, &energy);
  #endif
    exit(0);
  */

  /*
    if(lsms.global.checkIstop("buildKKRMatrix"))
    {
      if(lsms.global.iprint>=0)
      {
        FILE *f1=fopen("kkrmat.out","w");
        FILE *f2=fopen("kkrmat.pattern","w");
        for(int i=0; i<nrmat_ns; i++)
        {
          fprintf(f2,"%5d ",i);
          for(int j=0; j<nrmat_ns; j++)
          {
            fprintf(f1,"%5d %5d %lf
  %lf\n",i,j,std::real(m(i,j)),std::imag(m(i,j))); if(std::abs(m(i,j))<0.000001)
  fprintf(f2,"."); else fprintf(f2,"x");
          }
          fprintf(f2,"\n");
        }
        fclose(f1); fclose(f2);
      }
  // calculate matrix norm and condition number
      double work[2*nrmat_ns];
      int ipiv[nrmat_ns];
      int info;
      double norm=zlange_("1",&nrmat_ns,&nrmat_ns,&m(0,0),&nrmat_ns,work);
      printf("Matrix 1-norm of the kkr matrix=%lf\n",norm);
      zgetrf_(&nrmat_ns, &nrmat_ns, &m(0,0),&nrmat_ns,ipiv,&info);
      printf("ZGETRF on kkr-matrix returned info=%d\n",info);
      if(info==0)
      {
        zgetri_(&nrmat_ns,&m(0,0),&nrmat_ns,ipiv,(Complex
  *)work,&nrmat_ns,&info); printf("ZGETRI returned info=%d\n",info); double
  norm_i=zlange_("1",&nrmat_ns,&nrmat_ns,&m(0,0),&nrmat_ns,work); printf("Matrix
  1-norm of the kkr matrix inverse=%lf\n",norm_i); printf("1-norm condition
  number of the kkr-matrix=%lf\n",norm*norm_i);
      }
      // exit(0);
    }
  */

  delete[] gij;
  delete[] sinmp;
  delete[] cosmp;
  delete[] plm;
  delete[] hfn;
  delete[] dlm;
  delete[] bgij;
  delete[] tmat_n;
}

// calculateTauMatrix replaces gettaucl from LSMS_1. The communication is
// performed in calculateAllTauMatrices and the t matrices are in tmatStore,
// replacing vbig in gettaucl.
void calculateTauMatrix(LSMSSystemParameters &lsms, LocalTypeInfo &local,
                        AtomData &atom, int localAtomIndex, int ispin,
                        Complex energy, Complex prel, Complex *tau00_l,
                        Matrix<Complex> &m, int iie) {
#if MST_LINEAR_SOLVER_DEFAULT == 0x015 && !defined(USE_IRSXGESV)
  printf("Error: USE_IRSXGESV needs to be defined to use IRSXgesv\n");
  exit(1);
#endif
#if MST_LINEAR_SOLVER_DEFAULT == 0x014 && !defined(USE_XGETRF)
  printf("Error: USE_XGETRF needs to be defined to use IRSXgesv\n");
  exit(1);
#endif
  const unsigned int defaultLinearSolver = MST_LINEAR_SOLVER_DEFAULT;
  unsigned int linearSolver =
      lsms.global.linearSolver &
      MST_LINEAR_SOLVER_MASK;  // only use 12 least significant bits
  if (linearSolver == 0) linearSolver = defaultLinearSolver;

  unsigned int buildKKRMatrixKernel =
      lsms.global.linearSolver & MST_BUILD_KKR_MATRIX_MASK;
  if (buildKKRMatrixKernel == 0)
    buildKKRMatrixKernel = MST_BUILD_KKR_MATRIX_DEFAULT;

  int nrmat_ns = lsms.n_spin_cant * atom.nrmat;
  int kkrsz_ns = lsms.n_spin_cant * atom.kkrsz;

  Matrix<Complex> tau00(kkrsz_ns, kkrsz_ns);
  Complex *devM, *devT0;

  // =======================================
  // build the KKR matrix
  // =======================================
  m.resize(nrmat_ns, nrmat_ns);

  double timeBuildKKRMatrix = MPI_Wtime();

#ifdef USE_NVTX
  nvtxRangePushA("buildKKRMatrix");
#endif
  switch (buildKKRMatrixKernel) {
    case MST_BUILD_KKR_MATRIX_F77:
      buildKKRMatrix(lsms, local, atom, ispin, energy, prel, iie, m);
      break;
    case MST_BUILD_KKR_MATRIX_CPP:
      buildKKRMatrixCPU(lsms, local, atom, ispin, iie, energy, prel, m);
      break;
#if defined(ACCELERATOR_CUDA_C)
    case MST_BUILD_KKR_MATRIX_ACCELERATOR:
      /*
        {
      // test
      //  Matrix<Real> testLIZPos(3,atom.numLIZ);
      //  Matrix<Complex> bgij(nrmat_ns, nrmat_ns);
        Complex testIlp1[2*lsms.maxlmax + 1];
      //  cudaMemcpy(&bgij[0], devBgij, nrmat_ns*nrmat_ns*sizeof(Complex),
      cudaMemcpyDeviceToHost);
      //  cudaMemcpy(&testLIZPos[0], devAtom.LIZPos, 3*atom.numLIZ*sizeof(Real),
      cudaMemcpyDeviceToHost); cudaMemcpy(&testIlp1[0], DeviceConstants::ilp1,
      (2*lsms.maxlmax + 1)*sizeof(Complex), cudaMemcpyDeviceToHost);
        printf("calculateTauMatrix:\n");
        for(int l=0; l<2*lsms.maxlmax; l++)
        {
          printf("l=%d : ilp1 [%g + %gi] | DeviceConstats::ilp1 [%g +
      %gi]\n",l,IFactors::ilp1[l].real(),IFactors::ilp1[l].imag(),
      testIlp1[l].real(), testIlp1[l].imag());
        }
        }
      */
      devM = deviceStorage->getDevM();
      // printf("entering buildKKRMatrixCuda:\n");
      buildKKRMatrixCuda(lsms, local, atom, *deviceStorage,
                         deviceAtoms[localAtomIndex], ispin, iie, energy, prel,
                         devM);
      break;
#endif
#if defined(ACCELERATOR_HIP)
    case MST_BUILD_KKR_MATRIX_ACCELERATOR:
      /*
        {
      // test
      //  Matrix<Real> testLIZPos(3,atom.numLIZ);
      //  Matrix<Complex> bgij(nrmat_ns, nrmat_ns);
        Complex testIlp1[2*lsms.maxlmax + 1];
      //  deviceMemcpy(&bgij[0], devBgij, nrmat_ns*nrmat_ns*sizeof(Complex),
      deviceMemcpyDeviceToHost);
      //  deviceMemcpy(&testLIZPos[0], devAtom.LIZPos,
      3*atom.numLIZ*sizeof(Real), deviceMemcpyDeviceToHost);
        deviceMemcpy(&testIlp1[0], DeviceConstants::ilp1, (2*lsms.maxlmax +
      1)*sizeof(Complex), deviceMemcpyDeviceToHost);
        printf("calculateTauMatrix:\n");
        for(int l=0; l<2*lsms.maxlmax; l++)
        {
          printf("l=%d : ilp1 [%g + %gi] | DeviceConstats::ilp1 [%g +
      %gi]\n",l,IFactors::ilp1[l].real(),IFactors::ilp1[l].imag(),
      testIlp1[l].real(), testIlp1[l].imag());
        }
        }
      */
      devM = deviceStorage->getDevM();
      buildKKRMatrixHip(lsms, local, atom, *deviceStorage,
                        deviceAtoms[localAtomIndex], ispin, iie, energy, prel,
                        devM);
      break;
#endif
    default:
      printf("UNKNOWN KKR MARIX BUILD KERNEL (%x)!!!\n", buildKKRMatrixKernel);
      exit(1);
  }

  // -----------------------------------------------------------------------
  // make sure that the kkr matrix m is in the right place in memory
  // -----------------------------------------------------------------------
  if (lsms.lsmsMode == LSMSMode::kkrmat)
    linearSolver = MST_LINEAR_SOLVER_ZGETRF;  // make sure that the marix is on
                                              // the CPU side!

  switch (buildKKRMatrixKernel) {
    case MST_BUILD_KKR_MATRIX_F77:
    case MST_BUILD_KKR_MATRIX_CPP:
      // built on CPU:
      switch (linearSolver) {
        case MST_LINEAR_SOLVER_ZGESV:
        case MST_LINEAR_SOLVER_ZGETRF:
        case MST_LINEAR_SOLVER_ZCGESV:
        case MST_LINEAR_SOLVER_ZBLOCKLU_F77:
        case MST_LINEAR_SOLVER_ZBLOCKLU_CPP:
          break;
#if defined(ACCELERATOR_CUDA_C)
        case MST_LINEAR_SOLVER_ZGETRF_CUBLAS:
        case MST_LINEAR_SOLVER_ZBLOCKLU_CUBLAS:
        case MST_LINEAR_SOLVER_ZZGESV_CUSOLVER:
        case MST_LINEAR_SOLVER_ZGETRF_CUSOLVER:
        case MST_LINEAR_SOLVER_XGETRF_CUSOLVER:
        case MST_LINEAR_SOLVER_IRSXGESV_CUSOLVER:
          devM = deviceStorage->getDevM();
          transferMatrixToGPUCuda(devM, m);
          devT0 = deviceStorage->getDevT0();
          transferT0MatrixToGPUCuda(devT0, lsms, local, atom, iie, ispin);
          break;
#endif
#ifdef ACCELERATOR_HIP
        case MST_LINEAR_SOLVER_ZGETRF_ROCSOLVER:
          devM = deviceStorage->getDevM();
          transferMatrixToGPUHip(devM, m);
          devT0 = deviceStorage->getDevT0();
          transferT0MatrixToGPUHip(devT0, lsms, local, atom, iie);
          break;
#endif
        default:
          break;  // do nothing. We are using the CPU matrix
      }
      break;
#if defined(ACCELERATOR_CUDA_C)
    case MST_BUILD_KKR_MATRIX_ACCELERATOR:
      // built on GPU:
      switch (linearSolver) {
        case MST_LINEAR_SOLVER_ZGESV:
        case MST_LINEAR_SOLVER_ZGETRF:
        case MST_LINEAR_SOLVER_ZCGESV:
        case MST_LINEAR_SOLVER_ZBLOCKLU_F77:
        case MST_LINEAR_SOLVER_ZBLOCKLU_CPP:
          transferMatrixFromGPUCuda(m, (cuDoubleComplex *)devM);
          break;
        case MST_LINEAR_SOLVER_BLOCK_INVERSE_F77:
          transferMatrixFromGPUCuda(m, (cuDoubleComplex *)devM);
          break;
        case MST_LINEAR_SOLVER_BLOCK_INVERSE_CPP:
          transferMatrixFromGPUCuda(m, (cuDoubleComplex *)devM);
          break;
        case MST_LINEAR_SOLVER_ZGETRF_CUBLAS:
        case MST_LINEAR_SOLVER_ZBLOCKLU_CUBLAS:
        case MST_LINEAR_SOLVER_ZZGESV_CUSOLVER:
        case MST_LINEAR_SOLVER_ZGETRF_CUSOLVER:
        case MST_LINEAR_SOLVER_XGETRF_CUSOLVER:
        case MST_LINEAR_SOLVER_IRSXGESV_CUSOLVER:
          devT0 = deviceStorage->getDevT0();
          transferT0MatrixToGPUCuda(devT0, lsms, local, atom, iie, ispin);
          break;
#ifdef ACCELERATOR_HIP
        case MST_LINEAR_SOLVER_ZGETRF_ROCSOLVER:
          printf("MIXING HIP AND CUDA KERNELS (%x)!!!\n", buildKKRMatrixKernel);
          exit(1);
#endif
        default:
          break;  // do nothing. We are using the GPU matrix
      }
      break;
#endif
#if defined(ACCELERATOR_HIP)
    case MST_BUILD_KKR_MATRIX_ACCELERATOR:
      // built on GPU:
      switch (linearSolver) {
        case MST_LINEAR_SOLVER_ZGESV:
        case MST_LINEAR_SOLVER_ZGETRF:
        case MST_LINEAR_SOLVER_ZCGESV:
        case MST_LINEAR_SOLVER_ZBLOCKLU_F77:
        case MST_LINEAR_SOLVER_ZBLOCKLU_CPP:
          transferMatrixFromGPUHip(m, (deviceDoubleComplex *)devM);
          break;
        case MST_LINEAR_SOLVER_ZGETRF_ROCSOLVER:
          devT0 = deviceStorage->getDevT0();
          transferT0MatrixToGPUHip(devT0, lsms, local, atom, iie);
          break;
#ifdef ACCELERATOR_CUDA_C
        case MST_LINEAR_SOLVER_ZGETRF_CUBLAS:
        case MST_LINEAR_SOLVER_ZBLOCKLU_CUBLAS:
        case MST_LINEAR_SOLVER_ZZGESV_CUSOLVER:
        case MST_LINEAR_SOLVER_ZGETRF_CUSOLVER:
        case MST_LINEAR_SOLVER_XGETRF_CUSOLVER:
        case MST_LINEAR_SOLVER_IRSXGESV_CUSOLVER:
          printf("MIXING HIP AND CUDA KERNELS (%x)!!!\n", buildKKRMatrixKernel);
          exit(1);
          break;
#endif
        default:
          printf(
              "UNKNOWN SOLVER FOR MST_BUILD_KKR_MATRIX_ACCELERATOR (%x)!!!\n",
              linearSolver);
          exit(1);
          break;  // do nothing. We are using the GPU matrix
      }
      break;
#endif
    default:
      printf("UNKNOWN KKR MARIX BUILD KERNEL (%x)!!!\n", buildKKRMatrixKernel);
      exit(1);
  }
#ifdef USE_NVTX
  nvtxRangePop();
#endif

  timeBuildKKRMatrix = MPI_Wtime() - timeBuildKKRMatrix;
  if (lsms.global.iprint >= 0)
    printf("  timeBuildKKRMatrix=%lf\n", timeBuildKKRMatrix);

  // Write the kkrmatrix?
  if (lsms.lsmsMode == LSMSMode::kkrmat) {
    if (lsms.global.iprint >= 0) {
      char fname[256];
      snprintf(fname, 200, "kkrmat_%03d.dat",iie);
      FILE *f1 = fopen(fname, "w");
      snprintf(fname, 200, "kkrmat_%03d.h5",iie);
      hid_t h5_f1 = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

      snprintf(fname, 200, "tmat_%03d.dat",iie);
      FILE *f3 = fopen(fname, "w");
      snprintf(fname, 200, "tmat_%03d.h5",iie);
      hid_t h5_f3 = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

      Complex tmp;
      hid_t h5_complex = H5Tcreate(H5T_COMPOUND, sizeof(Complex));
      H5Tinsert (h5_complex, "real", 0, H5T_NATIVE_DOUBLE);
      H5Tinsert (h5_complex, "imag", sizeof(Real), H5T_NATIVE_DOUBLE);
      hsize_t kkrmat_dims[2]; kkrmat_dims[0] = kkrmat_dims[1] = nrmat_ns;
      hid_t h5_kkrmat_space = H5Screate_simple(2, kkrmat_dims, NULL);
      hid_t ds_kkrmat_r = H5Dcreate2(h5_f1, "kkrmat", h5_complex, h5_kkrmat_space,
                                     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); 

      H5Dwrite(ds_kkrmat_r, h5_complex, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m(0,0));
      H5Dclose(ds_kkrmat_r);
      H5Sclose(h5_kkrmat_space);

      hsize_t tmat_dims[2]; tmat_dims[0] = tmat_dims[1] = kkrsz_ns;
      hid_t h5_tmat_space = H5Screate_simple(2, tmat_dims, NULL);
      hid_t ds_tmat_r = H5Dcreate2(h5_f3, "tmat", h5_complex, h5_tmat_space,
                                     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      H5Dwrite(ds_tmat_r, h5_complex, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
         &local.tmatStore(iie * local.blkSizeTmatStore, atom.LIZStoreIdx[0]));
      H5Dclose(ds_tmat_r);
      H5Sclose(h5_tmat_space);

      // snprintf(fname, 200, "kkrmat_%03d.pattern",iie);
      // FILE *f2 = fopen(fname, "w");
      
      fprintf(f1, "# %6d x %6d\n", nrmat_ns, nrmat_ns);
      fprintf(f1, "# energy: %lg %lg\n", std::real(energy), std::imag(energy));
      fprintf(f3, "# %3d x %3d\n", kkrsz_ns, kkrsz_ns);
      fprintf(f3, "# energy: %lg %lg\n", std::real(energy), std::imag(energy));
      for (int i = 0; i < nrmat_ns; i++) {
        // fprintf(f2, "%5d ", i);
        for (int j = 0; j < nrmat_ns; j++) {
          fprintf(f1, "%6d %6d %lg %lg\n", i, j, std::real(m(i, j)),
                  std::imag(m(i, j)));
          // if (std::abs(m(i, j)) < 0.000001)
          //   fprintf(f2, ".");
          // else
          //   fprintf(f2, "x");
        }
        // fprintf(f2, "\n");
      }

      if (lsms.n_spin_pola == lsms.n_spin_cant) {  // non polarized or spin canted

        for (int i = 0; i < kkrsz_ns; i++) {
          for (int j = 0; j < kkrsz_ns; j++) {
            fprintf(f3, "%6d %6d %lg %lg\n", i, j, // tMatrix(i, j) =
              std::real(local.tmatStore(i + j * kkrsz_ns + iie * local.blkSizeTmatStore,
                            atom.LIZStoreIdx[0])),
              std::imag(local.tmatStore(i + j * kkrsz_ns + iie * local.blkSizeTmatStore,
                            atom.LIZStoreIdx[0])));
          }
        }

      } else {
        int jsm = kkrsz_ns * kkrsz_ns * ispin;

        for (int i = 0; i < kkrsz_ns; i++) {
          for (int j = 0; j < kkrsz_ns; j++) {
            fprintf(f3, "%6d %6d %lg %lg\n", i, j, // local.tmatStore(
              std::real(local.tmatStore(i + j * kkrsz_ns + iie * local.blkSizeTmatStore + jsm,
              atom.LIZStoreIdx[0])),
              std::imag(local.tmatStore(i + j * kkrsz_ns + iie * local.blkSizeTmatStore + jsm,
              atom.LIZStoreIdx[0])));
          }
        }
      }

      fclose(f1);
      H5Tclose(h5_complex);
      H5Fclose(h5_f1);
      // fclose(f2);
      fclose(f3);
      H5Fclose(h5_f3);
    }
    // exitLSMS(comm, 0);
  }

#ifdef USE_NVTX
  nvtxRangePushA("linearSolver");
#endif
  // use the new or old solvers?
  if (linearSolver <
      MST_LINEAR_SOLVER_BLOCK_INVERSE_F77)  // new solvers. Old solvers have
                                            // numbers > 0x8000. different
                                            // postpocessing required. 0 is the
                                            // default solver, for the time
                                            // being use the old LSMS_1.9 one
  {
    switch (linearSolver) {
      case MST_LINEAR_SOLVER_ZGESV:
        solveTau00zgesv(lsms, local, atom, iie, m, tau00, ispin);
        break;
      case MST_LINEAR_SOLVER_ZGETRF:
        solveTau00zgetrf(lsms, local, atom, iie, m, tau00, ispin);
        break;
#ifndef ARCH_IBM
      case MST_LINEAR_SOLVER_ZCGESV:
        solveTau00zcgesv(lsms, local, atom, iie, m, tau00, ispin);
        break;
#endif
      case MST_LINEAR_SOLVER_ZBLOCKLU_F77:
        solveTau00zblocklu_f77(lsms, local, atom, iie, m, tau00, ispin);
        break;
      case MST_LINEAR_SOLVER_ZBLOCKLU_CPP:
        solveTau00zblocklu_cpp(lsms, local, atom, iie, m, tau00, ispin);
        break;
#ifdef ACCELERATOR_CUDA_C
      case MST_LINEAR_SOLVER_ZGETRF_CUBLAS:
        solveTau00zgetrf_cublas(lsms, local, *deviceStorage, atom, devT0, devM,
                                tau00);
        break;
      case MST_LINEAR_SOLVER_ZBLOCKLU_CUBLAS:
        printf("MST_LINEAR_SOLVER_ZBLOCKLU_CUBLAS (%d) not implemented!!!\n",
               linearSolver);
        exit(1);
        break;
#ifndef ARCH_IBM
      case MST_LINEAR_SOLVER_ZZGESV_CUSOLVER:
        solveTau00zzgesv_cusolver(lsms, local, *deviceStorage, atom, devT0,
                                  devM, tau00, ispin);
        break;
#endif
      case MST_LINEAR_SOLVER_ZGETRF_CUSOLVER:
        solveTau00zgetrf_cusolver(lsms, local, *deviceStorage, atom, devT0,
                                  devM, tau00, ispin);
        break;
#ifdef USE_XGETRF
      case MST_LINEAR_SOLVER_XGETRF_CUSOLVER:
        solveTau00Xgetrf_cusolver(lsms, local, *deviceStorage, atom, devT0,
                                  devM, tau00, ispin);
        break;
#endif
#ifdef USE_IRSXGESV
      case MST_LINEAR_SOLVER_IRSXGESV_CUSOLVER:
        solveTau00IRSXgesv_cusolver(lsms, local, *deviceStorage, atom, devT0,
                                    devM, tau00, ispin);
        break;
#endif
#endif
#ifdef ACCELERATOR_HIP
      case MST_LINEAR_SOLVER_ZGETRF_ROCSOLVER:
        solveTau00zgetrf_rocsolver(lsms, local, *deviceStorage, atom, devT0,
                                   devM, tau00);
        break;
#endif
      default:
        printf("UNKNOWN LINEAR SOLVER (%d)!!!\n", linearSolver);
        exit(1);
    }

    double timePostproc = MPI_Wtime();
    if (lsms.relativity != full) {
      calculateTau00MinusT(lsms, local, atom, iie, tau00, tau00, ispin);
      rotateTau00ToLocalFrameNonRelativistic(lsms, atom, tau00, tau00_l);
    } else {
      rotateTau00ToLocalFrameRelativistic(lsms, atom, tau00, tau00_l);
    }
    timePostproc = MPI_Wtime() - timePostproc;
    if (lsms.global.iprint >= 1) printf("  timePostproc=%lf\n", timePostproc);

  } else {  // old solvers
    // if(!lsms.global.checkIstop("buildKKRMatrix"))
    {
      // invert matrix to get tau00
      // set up the block sizes for the block inversion:

      int nblk;
#if defined(ACCELERATOR_LIBSCI)
      nblk = 2;
#elif defined(ACCELERATOR_CUBLAS)
      // nblk=12;
      nblk = 4;
#elif defined(ACCELERATOR_CUDA_C)
      int max_blk_sz = 175;
      // assign blocks in a load balanced way
      if ((nrmat_ns - kkrsz_ns) % max_blk_sz == 0)
        nblk = (nrmat_ns - kkrsz_ns) / max_blk_sz + 1;
      else {
        nblk = (nrmat_ns - kkrsz_ns) / (max_blk_sz - 1) + 1;
        if ((nrmat_ns - kkrsz_ns) % (max_blk_sz - 1) > max_blk_sz / 2) nblk++;
      }
#else
      nblk = 4;
#endif
      if (kkrsz_ns == nrmat_ns)
        nblk = 1;
      else {
#if !defined(ACCELERATOR_CUDA_C) || defined(ACCELERATOR_CUBLAS)
        if (lsms.zblockLUSize > 0) {
          // assign blocks in a load balanced way
          if ((nrmat_ns - kkrsz_ns) % lsms.zblockLUSize == 0)
            nblk = (nrmat_ns - kkrsz_ns) / lsms.zblockLUSize + 1;
          else {
            nblk = (nrmat_ns - kkrsz_ns) / (lsms.zblockLUSize - 1) + 1;
            if ((nrmat_ns - kkrsz_ns) % (lsms.zblockLUSize - 1) >
                lsms.zblockLUSize / 2)
              nblk++;
          }
        }
#endif
      }
      int blk_sz[1000];
      assert(nblk <= 1000);

      blk_sz[0] = kkrsz_ns;
      if (nblk == 1)
        blk_sz[0] = nrmat_ns;
      else if (nblk == 2)
        blk_sz[1] = nrmat_ns - blk_sz[0];
      else if (nblk > 2)
      //  {
      //    int min_sz=(nrmat_ns-blk_sz[0])/(nblk-1);
      //    for(int i=1; i<nblk; i++) blk_sz[i]=min_sz;
      //    blk_sz[nblk-1]=nrmat_ns-blk_sz[0]-(nblk-2)*min_sz;
      //  }
      {
        int min_sz = (nrmat_ns - blk_sz[0]) / (nblk - 1);
        int rem = (nrmat_ns - blk_sz[0]) % (nblk - 1);
        int i = 1;
        for (; i <= rem; i++) blk_sz[i] = min_sz + 1;
        for (; i < nblk; i++) blk_sz[i] = min_sz;
      }
      if (lsms.global.iprint >= 1) {
        printf("nrmat_ns=%d\nnblk=%d\n", nrmat_ns, nblk);
        for (int i = 0; i < nblk; i++)
          printf("  blk_sz[%d]=%d\n", i, blk_sz[i]);
      }
      // block inversion:
      // vecs only needed for alg>2!
      // Complex vecs[nrmat_ns*(kkrsz_ns*6+6)];
      Complex *vecs = nullptr;

      int *ipvt = new int[nrmat_ns];
      Matrix<Complex> delta(kkrsz_ns, kkrsz_ns);
      // Complex *delta = new Complex[kkrsz_ns*kkrsz_ns];
      int *iwork = new int[kkrsz_ns * atom.numLIZ];
      Real *rwork = new Real[kkrsz_ns * atom.numLIZ];
      Complex *work1 = new Complex[kkrsz_ns * atom.numLIZ];
      int alg = 2;
      if (alg > 2)
        vecs = (Complex *)malloc((nrmat_ns * (kkrsz_ns * 6 + 6)) *
                                 sizeof(Complex));
      int *idcol = new int[blk_sz[0]];
      idcol[0] = 0;
      {
        if (linearSolver == MST_LINEAR_SOLVER_BLOCK_INVERSE_F77)
          block_inv_(&m(0, 0), vecs, &nrmat_ns, &nrmat_ns, &nrmat_ns, ipvt,
                     blk_sz, &nblk, &delta(0, 0), iwork, rwork, work1, &alg,
                     idcol, &lsms.global.iprint);
        else if (linearSolver == MST_LINEAR_SOLVER_BLOCK_INVERSE_CPP)
          block_inverse(m, blk_sz, nblk, delta, ipvt, idcol);
        else {
          printf("Unknown linear solver (%d)!!\n", linearSolver);
          exit(2);
        }
      }

      if (alg > 2) free(vecs);

      double timePostproc = MPI_Wtime();
      Matrix<Complex> tau00(kkrsz_ns, kkrsz_ns);
      if (lsms.relativity != full) {
        if (lsms.n_spin_pola ==
            lsms.n_spin_cant) {  // non polarized or spin canted

          tau_inv_postproc_nrel_(
              &kkrsz_ns, &lsms.n_spin_cant, &m(0, 0), &delta(0, 0),
              &local.tmatStore(iie * local.blkSizeTmatStore,
                               atom.LIZStoreIdx[0]),
              ipvt, &tau00(0, 0), atom.ubr, atom.ubrd, tau00_l);
        } else {
          int jsm = kkrsz_ns * kkrsz_ns * ispin;
          tau_inv_postproc_nrel_(
              &kkrsz_ns, &lsms.n_spin_cant, &m(0, 0), &delta(0, 0),
              &local.tmatStore(iie * local.blkSizeTmatStore + jsm,
                               atom.LIZStoreIdx[0]),
              ipvt, &tau00(0, 0), atom.ubr, atom.ubrd, tau00_l);
        }

      } else
        tau_inv_postproc_rel_(
            &kkrsz_ns, &m(0, 0), &delta(0, 0),
            &local.tmatStore(iie * local.blkSizeTmatStore, atom.LIZStoreIdx[0]),
            ipvt, &tau00(0, 0), &atom.dmat(0, 0), &atom.dmatp(0, 0), tau00_l);

      timePostproc = MPI_Wtime() - timePostproc;
      if (lsms.global.iprint >= 1) printf("  timePostproc=%lf\n", timePostproc);

      delete[] ipvt;
      // delete [] delta;
      delete[] iwork;
      delete[] rwork;
      delete[] work1;
      delete[] idcol;
    }
  }  // end of old solvers
#ifdef USE_NVTX
  nvtxRangePop();
#endif
}

// calculateAllTauMatrices replaces gettau and the communication part of
// gettaucl in LSMS_1
void calculateAllTauMatrices(
    LSMSCommunication &comm, LSMSSystemParameters &lsms, LocalTypeInfo &local,
    std::vector<Matrix<Real> > &vr, Complex energy, int iie,
    // std::vector<NonRelativisticSingleScattererSolution> &solution,
    Matrix<Complex> &tau00_l) {
  Complex prel = std::sqrt(energy * (1.0 + energy * c2inv));
  Complex pnrel = std::sqrt(energy);

  int max_nrmat_ns = 0;
  int max_kkrsz = 0;
  for (int i = 0; i < local.num_local; i++) {
    if (max_nrmat_ns < lsms.n_spin_cant * local.atom[i].nrmat)
      max_nrmat_ns = lsms.n_spin_cant * local.atom[i].nrmat;
    if (max_kkrsz < local.atom[i].kkrsz) max_kkrsz = local.atom[i].kkrsz;
  }

  /*
    std::vector<Matrix<Complex> >ms;

  #if !(defined _OPENMP) || (defined ACCELERATOR_CULA) ||
  defined(ACCELERATOR_LIBSCI) ms.resize(1); #if defined(ACCELERATOR_CULA) ||
  defined(ACCELERATOR_LIBSCI) int pinned; Complex *dat; pinned =
  cudaMallocHost((void**)&dat,max_nrmat_ns*max_nrmat_ns*sizeof(Complex)); if (
  pinned != cudaSuccess )
    {
      fprintf(stderr, "Matrix not pinned\n");
      dat = (Complex*)malloc(max_nrmat_ns*max_nrmat_ns*sizeof(Complex));
    }

    ms[0].retarget(max_nrmat_ns,max_nrmat_ns,dat);
  #else
    ms[0].resize(max_nrmat_ns,max_nrmat_ns);
  #endif
  #else
    ms.resize(omp_get_max_threads());
  #pragma omp parallel for default(none) shared(ms,max_nrmat_ns)
    for(int i=0; i<omp_get_max_threads(); i++)
      ms[i].resize(max_nrmat_ns,max_nrmat_ns);
  #endif
  */

  Complex *m_dat = NULL;
#if defined(ACCELERATOR_LIBSCI) || defined(ACCELERATOR_CUDA_C) || \
    defined(ACCELERATOR_HIP)
  m_dat = get_host_m_(max_nrmat_ns);
//  printf("m_dat = %p\n",m_dat);
#else
  m_dat = NULL;
#endif
  double timeCalcTauMatTotal = MPI_Wtime();
#if defined(ACCELERATOR_LIBSCI) || defined(ACCELERATOR_CUDA_C) || \
    defined(ACCELERATOR_HIP)
#pragma omp parallel for default(none) shared(                              \
    lsms, local, energy, prel, tau00_l, max_nrmat_ns, m_dat, deviceStorage) \
    firstprivate(iie) num_threads(lsms.global.GPUThreads)
#else
#pragma omp parallel for default(none) shared( \
    lsms, local, energy, prel, tau00_l, max_nrmat_ns, m_dat) firstprivate(iie)
#endif
  for (int i = 0; i < local.num_local; i++) {
    // printf("Num threads: %d\n",omp_get_num_threads());
    // printf("i: %d, threadId :%d\n",i,omp_get_thread_num());
#if defined(ACCELERATOR_LIBSCI) || defined(ACCELERATOR_CUDA_C) || \
    defined(ACCELERATOR_HIP)
    Matrix<Complex> m(
        max_nrmat_ns, max_nrmat_ns,
        m_dat + max_nrmat_ns * max_nrmat_ns * omp_get_thread_num());
#else
    Matrix<Complex> m(max_nrmat_ns, max_nrmat_ns);
#endif
    double timeCalcTauMat = MPI_Wtime();
    calculateTauMatrix(lsms, local, local.atom[i], i, 0, energy, prel,
                       &tau00_l(0, i), m, iie);
    if (lsms.n_spin_pola != lsms.n_spin_cant)  // spin polarized second spin
      calculateTauMatrix(lsms, local, local.atom[i], i, 1, energy, prel,
                         &tau00_l(0, i + local.num_local), m, iie);
    timeCalcTauMat = MPI_Wtime() - timeCalcTauMat;
    if (lsms.global.iprint >= 1)
      printf("calculateTauMatrix: %lf sec\n", timeCalcTauMat);
  }
  timeCalcTauMatTotal = MPI_Wtime() - timeCalcTauMatTotal;
  if (lsms.global.iprint >= 1)
    printf("calculateTauMatrix total time: %lf sec\n", timeCalcTauMatTotal);

/*
  ... various post processing
  and calculate green, dos, dosck...
*/
/*
  for(int i=0; i<ms.size(); i++)
    ms[i].resize(0,0);
#if defined(ACCELERATOR_CULA) || defined(ACCELERATOR_LIBSCI)
  if ( pinned == cudaSuccess )
    cudaFreeHost(dat);
  else
    free(dat);
#endif
*/
#if defined(ACCELERATOR_LIBSCI) || defined(ACCELERATOR_CUDA_C) || \
    defined(ACCELERATOR_HIP)

#else
  m_dat = NULL;
#endif
}
