/* -*- mode: C++; c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil
 * -*- */

#include <hip/hip_complex.h>
#include <hip/hip_runtime.h>
#include <hipblas/hipblas.h>
#include <rocsolver/rocsolver.h>
#include <stdio.h>

#include <vector>

#include "Accelerator/DeviceStorage.hpp"
#include "Complex.hpp"
#include "Matrix.hpp"
#include "linearSolvers.hpp"


/*
#define IDX(i, j, lDim) (((j)*(lDim))+(i))

template <typename T>
void zeroMatrixHip(T *devM, int lDim, int nCol)
{
//  for(int i=0; i<m.n_row(); i++)
//    for(int j=0; j<m.n_col(); j++)
//      m(i,j) = 0.0;
  hipError_t ret = hipMemset(devM, 0, lDim*nCol*sizeof(T));
}

template <typename T>
__global__ void setDiagonalKernelHip(T *devM, int lDim, int nCol, T val)
{
  int i=hipBlockIdx_x*hipBlockDim_x + hipThreadIdx_x;
  if(i<nCol)
  {
    devM[IDX(i, i, lDim)] = val;
  }
}

template <typename T>
__global__ void addDiagonalKernelHip(T *devM, int lDim, int nCol, T val)
{
  int i=hipBlockIdx_x*hipBlockDim_x + hipThreadIdx_x;
  if(i<nCol)
  {
    devM[IDX(i, i, lDim)] = hipCadd(devM[IDX(i, i, lDim)], val);
  }
}

template <typename T>
void unitMatrixHip(T *devM, int lDim, int nCol)
{
  hipError_t ret = zeroMatrixHip(devM, lDim, nCol);
  setDiagonalKernelHip<<<nCol,1>>>(devM, lDim, nCol, 1.0);
}
*/

template <typename T>
__global__ void zeroDiagonalBlocksKernelHip(T *devM, int lDim, int nCol,
                                            int blockSize) {
  int iBlock = hipBlockIdx_x * hipBlockDim_x + hipThreadIdx_x;
  int jBlock = hipBlockIdx_y * hipBlockDim_y + hipThreadIdx_y;
  if (iBlock < nCol / blockSize)
    if (jBlock < nCol / blockSize) {
      int ii = iBlock * blockSize;
      int jj = jBlock * blockSize;
      for (int i = 0; i < std::min(blockSize, nCol - ii); i++)
        for (int j = 0; j < std::min(blockSize, nCol - jj); j++)
          devM[IDX(ii + i, jj + j, lDim)] = 0.0;
    }
}

void transferT0MatrixToGPUHip(Complex *devT0, LSMSSystemParameters &lsms,
                              LocalTypeInfo &local, AtomData &atom, int iie) {
  int kkrsz_ns = lsms.n_spin_cant * atom.kkrsz;
  hipError_t ret = hipMemcpy(devT0,
            &local.tmatStore(iie * local.blkSizeTmatStore, atom.LIZStoreIdx[0]),
            kkrsz_ns * kkrsz_ns * sizeof(hipDoubleComplex),
            hipMemcpyHostToDevice);
}

void transferFullTMatrixToGPUHip(Complex *devT, LSMSSystemParameters &lsms, LocalTypeInfo &local,
                                  AtomData &atom, int ispin)
{
  int kkrsz_ns = lsms.n_spin_cant * atom.kkrsz;
  int nrmat_ns = lsms.n_spin_cant * atom.nrmat;
  Matrix <Complex> bigT(nrmat_ns, nrmat_ns);
  if (lsms.n_spin_pola == lsms.n_spin_cant) { // non polarized or spin canted
    for (int l = 0; l < atom.numLIZ; l++) {
      for (int i = 0; i < kkrsz_ns; i++) {
        for (int j = 0; j < kkrsz_ns; j++) {
          bigT(l*kkrsz_ns + i, l*kkrsz_ns + j) = local.tmatStore(i + j * kkrsz_ns, atom.LIZStoreIdx[l]);
        }
      }
    }
  } else {

    int jsm = kkrsz_ns * kkrsz_ns * ispin;
    for (int l = 0; l < atom.numLIZ; l++) {
      for (int i = 0; i < kkrsz_ns; i++) {
        for (int j = 0; j < kkrsz_ns; j++) {
          bigT(l*kkrsz_ns + i, l*kkrsz_ns + j) = local.tmatStore(i + j * kkrsz_ns + jsm, atom.LIZStoreIdx[l]);
        }
      }
    }
  }
  /*
  std::cout << "00 block of bigT on the CPU" << std::endl;
  for (int i = 0; i <kkrsz_ns; i++){
    for (int j = 0; j < kkrsz_ns; j++){
      std::cout << bigT(i,j) << "  ";
    }
    std::cout << std::endl;
  }
  */
  hipError_t ret = hipMemcpy(devT, &bigT(0,0), nrmat_ns*nrmat_ns*sizeof(hipDoubleComplex), hipMemcpyHostToDevice);
}

void transferMatrixToGPUHip(Complex *devM, Matrix<Complex> &m)
{
  hipError_t ret = hipMemcpy(devM, &m(0,0), m.l_dim()*m.n_col()*sizeof(hipDoubleComplex), hipMemcpyHostToDevice);
}

void transferMatrixFromGPUHip(Matrix<Complex> &m, hipDoubleComplex *devM) {
  hipError_t ret = hipMemcpy(&m(0, 0), devM, m.l_dim() * m.n_col() * sizeof(hipDoubleComplex),
            hipMemcpyDeviceToHost);
}

__global__ void copyTMatrixToTauHip(hipDoubleComplex *tau, hipDoubleComplex *t,
                                    int kkrsz, int nrmat) {
  int i = hipBlockIdx_x * hipBlockDim_x + hipThreadIdx_x;
  if (i < kkrsz) {
    for (int j = 0; j < kkrsz; j++) tau[IDX(i, j, nrmat)] = t[IDX(i, j, kkrsz)];
  }
}

__global__ void copyBigTMatrixToTauHip(hipDoubleComplex *tau, hipDoubleComplex *t, int nrmat)
{
   int i = hipBlockIdx_x*hipBlockDim_x + hipThreadIdx_x;
   if (i < nrmat)
   {
     for(int j=0;j<nrmat;j++){
       tau[IDX(i,j,nrmat)] = t[IDX(i,j,nrmat)];
     }
   }
}

__global__ void copyTauToTau00Hip(hipDoubleComplex *tau00, hipDoubleComplex *tau, int kkrsz, int nrmat)
{
  int i = hipBlockIdx_x*hipBlockDim_x + hipThreadIdx_x;
  if(i < kkrsz)
  {
    for(int j=0; j<kkrsz; j++)
      tau00[IDX(i,j,kkrsz)] = tau[IDX(i,j,nrmat)];
  }
}

/*
void solveTau00zgetrf_rocsolver(LSMSSystemParameters &lsms, LocalTypeInfo
&local, DeviceStorage &d, AtomData &atom, Complex *tMatrix, Complex *devM,
                             Matrix<Complex> &tau00)
{
  cublasHandle_t cublasHandle = DeviceStorage::getCublasHandle();
  int nrmat_ns = lsms.n_spin_cant*atom.nrmat; // total size of the kkr matrix
  int kkrsz_ns = lsms.n_spin_cant*atom.kkrsz; // size of t00 block
  // reference algorithm. Use LU factorization and linear solve for dense
matrices in LAPACK hipDoubleComplex *Aarray[1], *Barray[1];

  hipDoubleComplex *devTau = (hipDoubleComplex *)d.getDevTau();
  hipDoubleComplex *devTau00 = (hipDoubleComplex *)d.getDevTau00();
  // printf("zero Matrix\n");
  zeroMatrixHip(devTau, nrmat_ns, kkrsz_ns);
  // printf("copyTMatrixToTau\n");
  copyTMatrixToTauHip<<<kkrsz_ns,1>>>(devTau, (hipDoubleComplex *)tMatrix,
kkrsz_ns, nrmat_ns);

  Barray[0] = devTau;

  Aarray[0] = (hipDoubleComplex *)devM;

  int *ipivArray=d.getDevIpvt();
  int infoArray[1]; // d.getDevInfo();
  int info;

  // printf("cublasZgetrfBatched\n");
  cublasZgetrfBatched(cublasHandle, nrmat_ns, Aarray, nrmat_ns, ipivArray,
infoArray, 1);
  // printf("cublasZgetrsBatched\n");

  cublasZgetrsBatched(cublasHandle, CUBLAS_OP_N, nrmat_ns, kkrsz_ns, Aarray,
nrmat_ns, ipivArray, Barray, nrmat_ns, &info, 1);

  // copy result into tau00
  // printf("copyTauToTau00\n");
  copyTauToTau00Cuda<<<kkrsz_ns,1>>>(devTau00, devTau, kkrsz_ns, nrmat_ns);
  // printf("transferMatrixFromGPU\n");
  transferMatrixFromGPUCuda(tau00, devTau00);
}

#ifndef ARCH_IBM
void solveTau00zzgesv_cusolver(LSMSSystemParameters &lsms, LocalTypeInfo &local,
DeviceStorage &d, AtomData &atom, Complex *tMatrix, Complex *devM,
Matrix<Complex> &tau00)
{
  cusolverDnHandle_t cusolverDnHandle = DeviceStorage::getCusolverDnHandle();
  int nrmat_ns = lsms.n_spin_cant*atom.nrmat; // total size of the kkr matrix
  int kkrsz_ns = lsms.n_spin_cant*atom.kkrsz; // size of t00 block
  // reference algorithm. Use LU factorization and linear solve for dense
matrices in LAPACK

  cuDoubleComplex *devTau = (cuDoubleComplex *)d.getDevTau();
  cuDoubleComplex *devTau00 = (cuDoubleComplex *)d.getDevTau00();
  cuDoubleComplex *devWork = (cuDoubleComplex *)d.getDevWork();

  cuDoubleComplex *devT = (cuDoubleComplex *)d.getDevT();

  int *devIpiv = d.getDevIpvt();
  int devInfo[1]; // d.getDevInfo();

  zeroMatrixCuda(devTau, nrmat_ns, kkrsz_ns);
  zeroMatrixCuda(devT, nrmat_ns, kkrsz_ns);
  copyTMatrixToTauCuda<<<kkrsz_ns,1>>>(devT, (cuDoubleComplex *)tMatrix,
kkrsz_ns, nrmat_ns);

  int iter;

  cusolverStatus_t status = cusolverDnZZgesv(cusolverDnHandle, nrmat_ns,
kkrsz_ns, (cuDoubleComplex *)devM, nrmat_ns, devIpiv, devT, nrmat_ns, devTau,
nrmat_ns, devWork, d.getDevWorkBytes(), &iter, devInfo);

  if(status!=CUSOLVER_STATUS_SUCCESS)
  {
    printf("cusolverDnZZgesv returned %d\n",status);
  }

  copyTauToTau00Cuda<<<kkrsz_ns,1>>>(devTau00, devTau, kkrsz_ns, nrmat_ns);
  transferMatrixFromGPUCuda(tau00, devTau00);
}
#endif
*/

void solveTau00zgetrf_rocsolver(LSMSSystemParameters &lsms,
                                LocalTypeInfo &local, DeviceStorage &d,
                                AtomData &atom, Complex *tMatrix, Complex *devM,
                                Matrix<Complex> &tau00) {
  // cusolverDnHandle_t cusolverDnHandle = DeviceStorage::getCusolverDnHandle();
  rocblas_handle rocblasHandle = DeviceStorage::getRocBlasHandle();
  int nrmat_ns = lsms.n_spin_cant * atom.nrmat;  // total size of the kkr matrix
  int kkrsz_ns = lsms.n_spin_cant * atom.kkrsz;  // size of t00 block
  // reference algorithm. Use LU factorization and linear solve for dense
  // matrices in LAPACK
  hipDoubleComplex *devTau = (hipDoubleComplex *)d.getDevTau();
  hipDoubleComplex *devTau00 = (hipDoubleComplex *)d.getDevTau00();
  hipDoubleComplex *devWork = (hipDoubleComplex *)d.getDevWork();

  int *devIpiv = d.getDevIpvt();
  int *devInfo = d.getDevInfo();

  // printf("LSMS solveTau00zgetrf_rocsolver: entering\n");
  // fflush(stdout);
  zeroMatrixHip(devTau, nrmat_ns, kkrsz_ns);

#ifdef INCLUDE_ADDITIONAL_SYNCHRONIZE
  hipError_t ret = hipDeviceSynchronize();
#endif

  copyTMatrixToTauHip<<<kkrsz_ns, 1>>>(devTau, (hipDoubleComplex *)tMatrix,
                                       kkrsz_ns, nrmat_ns);

#ifdef INCLUDE_ADDITIONAL_SYNCHRONIZE
  ret = hipDeviceSynchronize();
#endif

  //  cusolverDnZgetrf(cusolverDnHandle, nrmat_ns, nrmat_ns,
  //                   (cuDoubleComplex *)devM, nrmat_ns, devWork, devIpiv,
  //                   devInfo );

  // printf("LSMS solveTau00zgetrf_rocsolver: before zgetrf\n");
  // fflush(stdout);
  rocsolver_zgetrf(rocblasHandle, nrmat_ns, nrmat_ns,
                   (rocblas_double_complex *)devM, nrmat_ns, devIpiv, devInfo);

#ifdef INCLUDE_ADDITIONAL_SYNCHRONIZE
  ret = hipDeviceSynchronize();
#endif

  //  cusolverDnZgetrs(cusolverDnHandle, CUBLAS_OP_N, nrmat_ns, kkrsz_ns,
  //                   (cuDoubleComplex *)devM, nrmat_ns, devIpiv, devTau,
  //                   nrmat_ns, devInfo);

  //  printf("LSMS solveTau00zgetrf_rocsolver: before zgetrs\n");
  // fflush(stdout);
  rocsolver_zgetrs(rocblasHandle, rocblas_operation_none, nrmat_ns, kkrsz_ns,
                   (rocblas_double_complex *)devM, nrmat_ns, devIpiv,
                   (rocblas_double_complex *)devTau, nrmat_ns);

#ifdef INCLUDE_ADDITIONAL_SYNCHRONIZE
  ret = hipDeviceSynchronize();
#endif

  // copy result into tau00
  // printf("LSMS solveTau00zgetrf_rocsolver: copy result into tau00\n");
  // fflush(stdout);
  copyTauToTau00Hip<<<kkrsz_ns, 1>>>(devTau00, devTau, kkrsz_ns, nrmat_ns);

#ifdef INCLUDE_ADDITIONAL_SYNCHRONIZE
  ret = hipDeviceSynchronize();
#endif

  transferMatrixFromGPUHip(tau00, devTau00);

#ifdef INCLUDE_ADDITIONAL_SYNCHRONIZE
  ret = hipDeviceSynchronize();
#endif

  // printf("LSMS solveTau00zgetrf_rocsolver: leaving\n");
  //  fflush(stdout);
}

void solveTauFullzgetrf_rocsolver(LSMSSystemParameters &lsms, LocalTypeInfo &local, DeviceStorage &d, AtomData &atom,
                               Complex *tMatrix, Complex *devM, Complex *devTauFull, int ispin)
{
  rocblas_handle rocblashandle = DeviceStorage::getRocBlasHandle();
  int nrmat_ns = lsms.n_spin_cant*atom.nrmat; // total size of the kkr matrix
  int kkrsz_ns = lsms.n_spin_cant*atom.kkrsz; // size of t00 block
  // reference algorithm. Use LU factorization and linear solve for dense matrices in LAPACK

  int *devIpiv = d.getDevIpvt();
  int *devInfo = d.getDevInfo();

  hipDoubleComplex *devWork = (hipDoubleComplex *)d.getDevWork();
  zeroMatrixHip(devTauFull, nrmat_ns, nrmat_ns);

#ifdef INCLUDE_ADDITIONAL_SYNCHRONIZE
  hipError_t ret = hipDeviceSynchronize();
#endif

  copyBigTMatrixToTauHip<<<nrmat_ns,1>>>((hipDoubleComplex *)devTauFull, (hipDoubleComplex *)tMatrix, nrmat_ns);

#ifdef INCLUDE_ADDITIONAL_SYNCHRONIZE
  ret = hipDeviceSynchronize();
#endif

  rocsolver_zgetrf(rocblashandle, nrmat_ns, nrmat_ns,
                                      (rocblas_double_complex *)devM, nrmat_ns, devIpiv, devInfo);

#ifdef INCLUDE_ADDITIONAL_SYNCHRONIZE
  ret = hipDeviceSynchronize();
#endif

  //std::cout << nrmat_ns << "  " << std::endl;
  //printf(" %p %p %p %p\n", devM, devTauFull, devIpiv, devInfo);
  rocsolver_zgetrs(rocblashandle, rocblas_operation_none, nrmat_ns, nrmat_ns,
                                      (rocblas_double_complex *)devM, nrmat_ns, devIpiv, (rocblas_double_complex *)devTauFull, nrmat_ns);

#ifdef INCLUDE_ADDITIONAL_SYNCHRONIZE
  ret = hipDeviceSynchronize();
#endif
  //transferMatrixFromGPUCuda(tau, devTauFull);
}
