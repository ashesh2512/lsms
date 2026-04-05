/* -*- mode: C++; c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil
 * -*- */

#include <stdio.h>

#include "Complex.hpp"
#include "Matrix.hpp"
#include <vector>

#include "Misc/Coeficients.hpp"
#include "Misc/Indices.hpp"

#include "Accelerator/DeviceStorage.hpp"
#include <hip/hip_runtime.h>
#include <hip/hip_complex.h>
#include <hip/hip_runtime.h>
#include <hipblas/hipblas.h>
#include <rocsolver/rocsolver.h>
#include <stdio.h>

#include <vector>

#include "Accelerator/DeviceStorage.hpp"
#include "Complex.hpp"
#include "Matrix.hpp"
#include "buildKKRMatrix.hpp"
#include "linearSolvers.hpp"

// we might want to distinguish between systems where all lmax (and consequently
// kkrsz_ns) are the same and systems with potential different lmax on different
// atoms and l steps

// #define COMPARE_ORIGINAL 1

// #ifdef COMPARE_ORIGINAL
// #define TOLERANCE 1e-06

// // if value is greater than 1, use relative value, else
// // use absolute value for comparison to avoid artificially big differences
// inline bool comp(float val, float test) {
//   bool comp = true;
//   if (std::abs(test) > 1.0) {
//     if (std ::abs((val - test) / test) > TOLERANCE) {
//       comp = false;
//     }
//   }
//   else {
//     if (std ::abs(val - test) > TOLERANCE) {
//       comp = false;
//     }
//   }
//   return comp;
// }
// #endif

// Fortran layout for matrix
// #define IDX(i, j, lDim) (((j)*(lDim))+(i))
#define IDX3(i, j, k, lDim, mDim) \
  (((k) * (lDim) * (mDim)) + ((j) * (lDim)) + (i))

__device__ __inline__ deviceFloatComplex complexExp_SP(deviceFloatComplex z) {
  float mExp = exp(hipCrealf(z));
  return (deviceFloatComplex) make_hipFloatComplex(mExp * cos(hipCimagf(z)),
                              mExp * sin(hipCimagf(z)));
}

__device__ inline void calculateHankelHip_SP(deviceFloatComplex prel, float r,
                                          int lend, deviceDoubleComplex *ilp1,
                                          deviceFloatComplex *hfn) {
  if (hipThreadIdx_x == 0) {
    const deviceFloatComplex sqrtm1 =  make_hipFloatComplex(0.0, 1.0);
    deviceFloatComplex z = hipCmulf(prel , make_hipFloatComplex(r, 0.0) );

    hfn[0] = make_hipFloatComplex(0.0, -1.0);
    hfn[1] = make_hipFloatComplex(-1.0, 0.0) - hipCdivf(sqrtm1, z);
    for (int l = 1; l < lend; l++) {
      hfn[l + 1] = hipCdivf((2.0 * l + 1.0) * hfn[l], z) - hfn[l - 1];
    }

    z = complexExp_SP( hipCmulf(sqrtm1, z)) / r;
    for (int l = 0; l <= lend; l++) {
      deviceFloatComplex ilp1_f    = make_hipFloatComplex((float)hipCreal(ilp1[l]), (float)hipCimag(ilp1[l]));
      hfn[l] = hipCmulf( hipCmulf((-hfn[l]), z), ilp1_f );
    }
  }
  //  __syncthreads();
}

__device__ inline void calculateSinCosPowersHip_SP(float *rij, int lend,
                                                float *sinmp, float *cosmp) {
  const float ptol = 1.0e-6;
  float pmag = sqrt(rij[0] * rij[0] + rij[1] * rij[1]);
  cosmp[0] = 1.0;
  sinmp[0] = 0.0;
  if (pmag > ptol) {
    cosmp[1] = rij[0] / pmag;
    sinmp[1] = rij[1] / pmag;
  } else {
    cosmp[1] = 0.0;
    sinmp[1] = 0.0;
  }
  for (int m = 2; m <= lend; m++) {
    cosmp[m] = cosmp[m - 1] * cosmp[1] - sinmp[m - 1] * sinmp[1];
    sinmp[m] = sinmp[m - 1] * cosmp[1] + cosmp[m - 1] * sinmp[1];
  }
}

#define PLM_IDX(l, m) (((l) * ((l) + 1)) / 2 + (m))

__device__ __inline__ int plmIdxDev(int l, int m) {
  return l * (l + 1) / 2 + m;
}

__device__ void associatedLegendreFunctionNormalizedHip_SP(float x, int lmax,
                                                        float *Plm) {
  const float pi = acos(-1.0);
  // y = \sqrt{1-x^2}
  float y = sqrt(1.0 - x * x);
  // initialize the first entry
  Plm[0] = sqrt(1.0 / (4.0 * pi));

  if (lmax < 1) return;

  for (int m = 1; m <= lmax; m++) {
    // \bar{P}_{mm} = - \sqrt{\frac{2m+1}{2m}} y \bar{P}_{m-1, m-1}
    Plm[plmIdxDev(m, m)] =
        -sqrt(float(2 * m + 1) / float(2 * m)) * y * Plm[plmIdxDev(m - 1, m - 1)];
    // \bar{P}_{mm-1} = \sqrt{2 m + 1} x \bar{P}_{m-1, m-1}
    Plm[plmIdxDev(m, m - 1)] =
        sqrt(float(2 * m + 1)) * x * Plm[plmIdxDev(m - 1, m - 1)];
  }

  for (int m = 0; m < lmax; m++) {
    for (int l = m + 2; l <= lmax; l++) {
      // \bar{P}_{lm} = a_{lm} (x \bar{P}_{l-1. m} - b_{lm} \bar{P}_{l-2, m})
      // a_{lm} = \sqrt{\frac{(4 l^2 - 1)(l^2 - m^2)}}
      // b_{lm} = \sqrt{\frac{(l -1)^2 - m^2}{4 (l-1)^2 -1}}
      float a_lm = sqrt(float(4 * l * l - 1) / float(l * l - m * m));
      float b_lm = sqrt(float((l - 1) * (l - 1) - m * m) /
                        float(4 * (l - 1) * (l - 1) - 1));
      Plm[plmIdxDev(l, m)] = a_lm * (x * Plm[plmIdxDev(l - 1, m)] -
                                     b_lm * Plm[plmIdxDev(l - 2, m)]);
    }
  }
}

__device__ __inline__ deviceFloatComplex dlmFunction_SP(deviceFloatComplex *hfn,
                                                     float *cosmp,
                                                     float *sinmp,
                                                     float *plm, int l,
                                                     int m) {
  int mAbs = abs(m);

  deviceFloatComplex dlm = hfn[l] * plm[plmIdxDev(l, mAbs)];
  if (m == 0) return dlm;

  if (m < 0)
  {
    dlm = hipCmulf(dlm, make_hipFloatComplex(cosmp[mAbs], sinmp[mAbs]));

    if ((mAbs & 0x01) != 0) // m is odd
      dlm = -dlm;
  }
  else
  {
    dlm = hipCmulf(dlm, make_hipFloatComplex(cosmp[mAbs], -sinmp[mAbs]));
  }

  return dlm;
}

size_t sharedMemoryBGijHip_SP(LSMSSystemParameters &lsms, size_t *hfnOffset,
                           size_t *sinmpOffset, size_t *cosmpOffset,
                           size_t *plmOffset, size_t *dlmOffset) {
  size_t size = 0;

  *hfnOffset = size;
  size += sizeof(deviceFloatComplex) * (2 * lsms.maxlmax + 1);

  *sinmpOffset = size;
  size += sizeof(float) * (2 * lsms.maxlmax + 1);

  *cosmpOffset = size;
  size += sizeof(float) * (2 * lsms.maxlmax + 1);

  *plmOffset = size;
  size += sizeof(float) * (AngularMomentumIndices::ndlm);

  return size;
}

__global__ void setBGijHip_SP(bool fullRelativity, int n_spin_cant, int *LIZlmax,
                           int *offsets, size_t nrmat_ns,
                           deviceFloatComplex *devBgij) {
  if (n_spin_cant == 1) return;

  int ir1 = hipBlockIdx_x;
  int ir2 = hipBlockIdx_y;
  int iOffset = offsets[ir1];
  int jOffset = offsets[ir2];

  int kkri = (LIZlmax[ir1] + 1) * (LIZlmax[ir1] + 1);
  int kkrj = (LIZlmax[ir2] + 1) * (LIZlmax[ir2] + 1);

  if (!fullRelativity)  //(lsms.relativity != full)
  {
    for (int ij = hipThreadIdx_x; ij < kkri * kkrj; ij += hipBlockDim_x) {
      int i = ij % kkri;
      int j = ij / kkri;

      devBgij[IDX(iOffset + kkri + i, jOffset + j, nrmat_ns)] =
          make_hipFloatComplex(0.0, 0.0);  // bgij(iOffset + i, jOffset + j);
      devBgij[IDX(iOffset + i, jOffset + kkrj + j, nrmat_ns)] =
          make_hipFloatComplex(0.0, 0.0);  // bgij(iOffset + i, jOqffset + j);
      devBgij[IDX(iOffset + kkri + i, jOffset + kkrj + j, nrmat_ns)] =
          devBgij[IDX(iOffset + i, jOffset + j, nrmat_ns)];
    }
  } else {
    printf(
        "Fully relativistic calculation not yet implemented in "
        "'MultipleScattering/buildKKRMatrix.cpp : setBGijCPU'\n");
    // exit(1);
  }
}

__global__ void buildGijHipKernel_SP(
    Real *LIZPos, int *LIZlmax, int *lofk, int *mofk, deviceDoubleComplex *ilp1,
    deviceDoubleComplex *illp, Real *cgnt, int ndlj_illp, int lmaxp1_cgnt,
    int ndlj_cgnt, size_t hfnOffset, size_t sinmpOffset, size_t cosmpOffset,
    size_t plmOffset, size_t dlmOffset,
// #if !defined(COMPARE_ORIGINAL)
    deviceFloatComplex energy, deviceFloatComplex prel, int *offsets,
    size_t nrmat_ns, deviceFloatComplex *devBgij)
// #else
//     deviceFloatComplex energy, deviceFloatComplex prel, int *offsets,
//     size_t nrmat_ns, deviceFloatComplex *devBgij, char *testSM)
// #endif
{
  int ir1 = hipBlockIdx_x;
  int ir2 = hipBlockIdx_y;

  HIP_DYNAMIC_SHARED(char, sharedMemory);

  if (ir1 != ir2) {
    int iOffset = offsets[ir1];
    // int iOffset = ir1 * kkrsz_ns;
    int jOffset = offsets[ir2];
    // int jOffset = ir2 * kkrsz_ns;

    float rij[3];
    rij[0] = LIZPos[3 * ir1 + 0] - LIZPos[3 * ir2 + 0];
    rij[1] = LIZPos[3 * ir1 + 1] - LIZPos[3 * ir2 + 1];
    rij[2] = LIZPos[3 * ir1 + 2] - LIZPos[3 * ir2 + 2];

    deviceFloatComplex *hfn =
        (deviceFloatComplex *)(sharedMemory + hfnOffset);
    float *sinmp = (float *)(sharedMemory + sinmpOffset);
    float *cosmp = (float *)(sharedMemory + cosmpOffset);
    float *plm = (float *)(sharedMemory + plmOffset);

// #if defined(COMPARE_ORIGINAL)
//     deviceDoubleComplex *testHfn = (deviceDoubleComplex *)(testSM + hfnOffset);
//     Real *testSinmp = (Real *)(testSM + sinmpOffset);
//     Real *testCosmp = (Real *)(testSM + cosmpOffset);
//     Real *testPlm = (Real *)(testSM + plmOffset);
//     // deviceDoubleComplex *testDlm = (deviceDoubleComplex *) (testSM +
//     // dlmOffset);
// #endif

    float r = sqrt(rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2]);
    int lmax1 = LIZlmax[ir1];
    int lmax2 = LIZlmax[ir2];
    int kkri = (lmax1 + 1) * (lmax1 + 1);
    int kkrj = (lmax2 + 1) * (lmax2 + 1);
    int lend = lmax1 + lmax2;

    const float pi4 = 4.0 * 2.0 * asin(1.0);
    float cosTheta = rij[2] / r;

    if (hipThreadIdx_x == 0) {
      calculateHankelHip_SP(prel, r, lend, ilp1, hfn);

      associatedLegendreFunctionNormalizedHip_SP(cosTheta, lend, plm);
 
      calculateSinCosPowersHip_SP(rij, lend, sinmp, cosmp);
    }
    __syncthreads();


// #if defined(COMPARE_ORIGINAL)
//     if (ir1 == 0 && ir2 == 1 && hipThreadIdx_x == 0) {
//       for (int l = 0; l <= lend; l++) {
//         testHfn[l] = hfn[l];
//         testSinmp[l] = sinmp[l];
//         testCosmp[l] = cosmp[l];
//       }
//     }
// #endif

    //     ================================================================
    //     calculate g(R_ij)...............................................
    for (int ij = hipThreadIdx_x; ij < kkri * kkrj; ij += hipBlockDim_x) {
      int lm2 = ij % kkri;
      int lm1 = ij / kkri;

      auto devBgijValue = make_hipFloatComplex(0.0, 0.0);

      int l1 = lofk[lm1];
      int m1 = mofk[lm1];

      int l2 = lofk[lm2];
      int m2 = mofk[lm2];

      int m3 = m2 - m1;
      int llow = max(abs(m3), abs(l1 - l2));
      if (hipCabsf(prel) == 0.0) llow = l1 + l2;

      int cgntIdx0 = IDX3(0, lm1, lm2, lmaxp1_cgnt, ndlj_cgnt);

      for (int l3 = l1 + l2; l3 >= llow; l3 -= 2) {
        int j = l3 * (l3 + 1) + m3;
   
        devBgijValue =
            devBgijValue
            + (float)cgnt[(l3 >> 1) + cgntIdx0] *
                  dlmFunction_SP(hfn, cosmp, sinmp, plm, l3, m3);  // dlm[j];
      }

      deviceFloatComplex illp_f = 
        make_hipFloatComplex((float)hipCreal(illp[IDX(lm2, lm1, ndlj_illp)]), (float)hipCimag(illp[IDX(lm2, lm1, ndlj_illp)]));

      devBgij[IDX(iOffset + lm2, jOffset + lm1, nrmat_ns)] =
        pi4 * hipCmulf(devBgijValue, illp_f);
    }
  }
}

__global__ void buildKKRMatrixMultiplyKernelHip_SP(
    int *LIZlmax, int *LIZStoreIdx, int *offsets, int kkrsz_ns, int ispin,
    int n_spin_pola, int n_spin_cant, int iie, int blkSizeTmatStore,
    int tmatStoreLDim, deviceFloatComplex *devTmatStore, int nrmat_ns,
    deviceFloatComplex *devBgij, deviceFloatComplex *devM) {
  int ir1 = hipBlockIdx_x;
  int ir2 = hipBlockIdx_y;

  deviceFloatComplex *tmat_n;
  int iOffset = offsets[ir1];
  int jOffset = offsets[ir2];

  if (ir1 != ir2) {
    int lmax1 = LIZlmax[ir1];
    int lmax2 = LIZlmax[ir2];
    int kkr1 = (lmax1 + 1) * (lmax1 + 1);
    int kkr2 = (lmax2 + 1) * (lmax2 + 1);
    int kkr1_ns = kkr1 * n_spin_cant;
    int kkr2_ns = kkr2 * n_spin_cant;

    tmat_n = &devTmatStore[IDX(iie * blkSizeTmatStore, LIZStoreIdx[ir1],
                               tmatStoreLDim)];

    for (int ij = hipThreadIdx_x; ij < kkr1_ns * kkr2_ns; ij += hipBlockDim_x) {
      int i = ij % kkr1_ns;
      int j = ij / kkr1_ns;

      auto devMValue = make_hipFloatComplex(0.0, 0.0);
      int devBgijIdx0 = IDX(iOffset, jOffset + j, nrmat_ns);
      int tmat_nIdx = IDX(i, 0, kkr1_ns);
      for (int k = 0; k < kkr1_ns; k++)
      {
        devMValue = devMValue -
                    hipCmulf(tmat_n[tmat_nIdx], devBgij[k + devBgijIdx0]);
        tmat_nIdx += kkr1_ns;
      }
      devM[IDX(iOffset + i, jOffset + j, nrmat_ns)] = devMValue;
    }
  }
}

__global__ void buildKKRMatrixMultiplyKernelHip_collinear_SP(
    int *LIZlmax, int *LIZStoreIdx, int *offsets, int kkrsz, int ispin,
    int iie, int blkSizeTmatStore,
    int tmatStoreLDim, deviceFloatComplex *devTmatStore, int nrmat_ns,
    deviceFloatComplex *devBgij, deviceFloatComplex *devM) {
  int ir1 = hipBlockIdx_x;
  int ir2 = hipBlockIdx_y;

  deviceFloatComplex *tmat_n;
  int iOffset = offsets[ir1];
  int jOffset = offsets[ir2];
  int spinOffset = kkrsz * kkrsz * ispin;

  if (ir1 != ir2) {
    int lmax1 = LIZlmax[ir1];
    int lmax2 = LIZlmax[ir2];
    int kkr1 = (lmax1 + 1) * (lmax1 + 1);
    int kkr2 = (lmax2 + 1) * (lmax2 + 1);
    int kkr1_ns = kkr1;
    int kkr2_ns = kkr2;

    tmat_n = &devTmatStore[IDX(iie * blkSizeTmatStore + spinOffset, LIZStoreIdx[ir1],
                               tmatStoreLDim)];

    for (int ij = hipThreadIdx_x; ij < kkr1_ns * kkr2_ns; ij += hipBlockDim_x) {
      int i = ij % kkr1_ns;
      int j = ij / kkr1_ns;

      auto devMValue = make_hipFloatComplex(0.0, 0.0);
      int devBgijIdx0 = IDX(iOffset, jOffset + j, nrmat_ns);
      int tmat_nIdx = IDX(i, 0, kkr1_ns);
      for (int k = 0; k < kkr1_ns; k++)
      {

        devMValue = devMValue -
                    hipCmulf(tmat_n[tmat_nIdx], devBgij[k + devBgijIdx0]);

        tmat_nIdx += kkr1_ns;
      }
      devM[IDX(iOffset + i, jOffset + j, nrmat_ns)] = devMValue;
    }
  }
}

void buildKKRMatrixLMaxIdenticalHip_SP(LSMSSystemParameters &lsms,
                                    LocalTypeInfo &local, AtomData &atom,
                                    DeviceStorage &d, DeviceAtom &devAtom,
                                    int ispin, int iie, Complex energy,
                                    Complex prel, ComplexF *devM) {
  hipblasHandle_t hipblasHandle = DeviceStorage::getHipBlasHandle();
  int nrmat_ns = lsms.n_spin_cant * atom.nrmat;  // total size of the kkr matrix
  int kkrsz_ns = lsms.n_spin_cant * atom.kkrsz;  // size of t00 block
  bool fullRelativity = false;
  if (lsms.relativity == full) fullRelativity = true;

  ComplexF *devBgij = d.getDevBGijF();

  deviceFloatComplex cuEnergy = make_hipFloatComplex(energy.real(), energy.imag());
  deviceFloatComplex cuPrel = make_hipFloatComplex(prel.real(), prel.imag());

  unitMatrixHip<ComplexF>(devM, nrmat_ns, nrmat_ns);
  zeroMatrixHip(devBgij, nrmat_ns, nrmat_ns);

  // calculate Bgij
  // reuse ipvt for offsets
  int *devOffsets = d.getDevIpvt();

  std::vector<int> offsets(devAtom.numLIZ);
  for (int ir = 0; ir < devAtom.numLIZ; ir++) offsets[ir] = ir * kkrsz_ns;

  deviceError_t ret = deviceMemcpy(devOffsets, &offsets[0], atom.numLIZ * sizeof(int),
               deviceMemcpyHostToDevice);

  size_t hfnOffset, sinmpOffset, cosmpOffset, plmOffset, dlmOffset;
  size_t smSize = sharedMemoryBGijHip_SP(lsms, &hfnOffset, &sinmpOffset,
                                      &cosmpOffset, &plmOffset, &dlmOffset);
// #ifdef COMPARE_ORIGINAL
//   char *devTestSM;
//   ret = deviceMalloc((void **)&devTestSM, smSize);
// #endif
  // int threads = 256;
  int threads = 1;
  dim3 blocks = dim3(devAtom.numLIZ, devAtom.numLIZ, 1);
  buildGijHipKernel_SP<<<blocks, threads, smSize>>>(
      devAtom.LIZPos, devAtom.LIZlmax, DeviceConstants::lofk,
      DeviceConstants::mofk, DeviceConstants::ilp1, DeviceConstants::illp,
      DeviceConstants::cgnt, DeviceConstants::ndlj_illp,
      DeviceConstants::lmaxp1_cgnt, DeviceConstants::ndlj_cgnt, hfnOffset,
      sinmpOffset, cosmpOffset, plmOffset, dlmOffset, cuEnergy, cuPrel,
// #if !defined(COMPARE_ORIGINAL)
      devOffsets, nrmat_ns, (deviceFloatComplex *)devBgij);
// #else
//       devOffsets, nrmat_ns, (deviceFloatComplex *)devBgij, devTestSM);
// #endif

  setBGijHip_SP<<<blocks, threads>>>(fullRelativity, lsms.n_spin_cant,
                                  devAtom.LIZlmax, devOffsets, nrmat_ns,
                                  (deviceFloatComplex *)devBgij);

// #ifdef COMPARE_ORIGINAL
//   Matrix<Real> testLIZPos(3, atom.numLIZ);
//   Matrix<Complex> bgij(nrmat_ns, nrmat_ns);
//   Complex testIlp1[2 * lsms.maxlmax + 1];
//   ret = deviceMemcpy(&bgij[0], devBgij, nrmat_ns * nrmat_ns * sizeof(Complex),
//                deviceMemcpyDeviceToHost);
//   ret = deviceMemcpy(&testLIZPos[0], devAtom.LIZPos, 3 * atom.numLIZ * sizeof(Real),
//                deviceMemcpyDeviceToHost);
//   ret = deviceMemcpy(&testIlp1[0], DeviceConstants::ilp1,
//                (2 * lsms.maxlmax + 1) * sizeof(Complex),
//                deviceMemcpyDeviceToHost);

//   for (int l = 0; l < 2 * lsms.maxlmax; l++) {
//     printf("l=%d : ilp1 [%g + %gi] | DeviceConstats::ilp1 [%g + %gi]\n", l,
//            IFactors::ilp1[l].real(), IFactors::ilp1[l].imag(),
//            testIlp1[l].real(), testIlp1[l].imag());
//   }

//   Complex testHfn[2 * lsms.maxlmax + 1];
//   Real testSinmp[2 * lsms.maxlmax + 1];
//   Real testCosmp[2 * lsms.maxlmax + 1];
//   // Real plm[((lsms.maxlmax+1) * (lsms.maxlmax+2)) / 2];
//   Real testPlm[AngularMomentumIndices::ndlm];
// // Complex testDlm[AngularMomentumIndices::ndlj];
//   ret = deviceMemcpy(testHfn, devTestSM + hfnOffset, (2*lsms.maxlmax + 1)*sizeof(Complex), deviceMemcpyDeviceToHost);
//   ret = deviceMemcpy(testSinmp, devTestSM + sinmpOffset, (2*lsms.maxlmax + 1)*sizeof(Real), deviceMemcpyDeviceToHost);
//   ret = deviceMemcpy(testCosmp, devTestSM + cosmpOffset, (2*lsms.maxlmax + 1)*sizeof(Real), deviceMemcpyDeviceToHost);
//   ret = deviceMemcpy(testPlm, devTestSM + plmOffset, AngularMomentumIndices::ndlm*sizeof(Real), deviceMemcpyDeviceToHost);
// // ret = deviceMemcpy(testDlm, devTestSM + dlmOffset, AngularMomentumIndices::ndlj*sizeof(Complex), deviceMemcpyDeviceToHost);

//   for (int i = 0; i < atom.numLIZ; i++) {
//     if (atom.LIZPos(0, i) != testLIZPos(0, i) ||
//         atom.LIZPos(1, i) != testLIZPos(1, i) ||
//         atom.LIZPos(2, i) != testLIZPos(2, i)) {
//       printf("atom.LIZPos(*,%d) [%.15f,%.15f,%.15f] != devAtom.LIZPos(*,%d) [%.15f,%.15f,%.15f]\n",
//         i,atom.LIZPos(0,i),atom.LIZPos(1,i),atom.LIZPos(2,i),
//         i,testLIZPos(0,i),testLIZPos(1,i),testLIZPos(2,i));
//     }
//   }
//   // loop over the LIZ blocks
//   Complex hfn[2 * lsms.maxlmax + 1];
//   Real sinmp[2 * lsms.maxlmax + 1];
//   Real cosmp[2 * lsms.maxlmax + 1];
//   // Real plm[((lsms.maxlmax+1) * (lsms.maxlmax+2)) / 2];
//   Real plm[AngularMomentumIndices::ndlm];
//   Complex dlm[AngularMomentumIndices::ndlj];
//   Real rij[3];
//   Real pi4 = 4.0 * 2.0 * std::asin(1.0);
//   bool exitCompare = false;
//   for (int ir1 = 0; ir1 < atom.numLIZ; ir1++) {
//     int iOffset =
//         ir1 *
//         kkrsz_ns;  // this assumes that there are NO lStep reductions of lmax!!!
//     for (int ir2 = 0; ir2 < atom.numLIZ; ir2++) {
//       int jOffset = ir2 * kkrsz_ns;  // this assumes that there are NO lStep
//                                      // reductions of lmax!!
//       int lmax1 = atom.LIZlmax[ir1];
//       int lmax2 = atom.LIZlmax[ir2];
//       int kkri = (lmax1 + 1) * (lmax1 + 1);
//       int kkrj = (lmax2 + 1) * (lmax2 + 1);
//       rij[0] = atom.LIZPos(0, ir1) - atom.LIZPos(0, ir2);
//       rij[1] = atom.LIZPos(1, ir1) - atom.LIZPos(1, ir2);
//       rij[2] = atom.LIZPos(2, ir1) - atom.LIZPos(2, ir2);
//       if (ir1 != ir2) {
//         int kkr1 = kkri;
//         int kkr2 = kkrj;
//         // bool exitCompare = false;
//         Matrix<Complex> gijTest(kkr1, kkr2);
//         Matrix<Complex> bgijTest(2 * kkr1, 2 * kkr2);
//         int lmax = lsms.maxlmax;
//         int kkrsz = (lmax + 1) * (lmax + 1);
//         makegij_(&atom.LIZlmax[ir1], &kkr1, &atom.LIZlmax[ir2], &kkr2,
//                  &lsms.maxlmax, &kkrsz, &AngularMomentumIndices::ndlj,
//                  &AngularMomentumIndices::ndlm, &prel, &rij[0], &sinmp[0],
//                  &cosmp[0], &SphericalHarmonicsCoeficients::clm[0], &plm[0],
//                  &GauntCoeficients::cgnt(0, 0, 0), &GauntCoeficients::lmax,
//                  &AngularMomentumIndices::lofk[0],
//                  &AngularMomentumIndices::mofk[0], &IFactors::ilp1[0],
//                  &IFactors::illp(0, 0), &hfn[0], &dlm[0], &gijTest(0, 0), &pi4,
//                  &lsms.global.iprint, lsms.global.istop, 32);

//         if (ir1 == 0 && ir2 == 1) {
//           for (int l = 0; l <= atom.LIZlmax[ir1] + atom.LIZlmax[ir2]; l++) {
//             if (!comp(sinmp[l], testSinmp[l]))
//               printf("sinmp[%d] (%.15f) != testSinmp[%d] (%.15f)\n", l, sinmp[l], l, testSinmp[l]);
//             if (!comp(cosmp[l], testCosmp[l]))
//               printf("cosmp[%d] (%.15f) != testCosmp[%d] (%.15f)\n", l, cosmp[l], l, testCosmp[l]);
//             if (!comp(hfn[l].real(), testHfn[l].real()) ||
//                 !comp(hfn[l].imag(), testHfn[l].imag()))
//               printf("hfn[%d] (%.15f + %.15fi) != testHfn[%d] (%.15f + %.15fi)\n", l, hfn[l].real(), hfn[l].imag(), l, testHfn[l].real(), testHfn[l].imag());
//           }
//         }

//         int idx = 0;
//         for (int i = 0; i < kkri; i++)
//           for (int j = 0; j < kkrj; j++) {
//             if (!comp(bgij(iOffset + i, jOffset + j).real(), gijTest(i, j).real()) ||
//                 !comp(bgij(iOffset + i, jOffset + j).imag(), gijTest(i, j).imag()))
//             {
//               printf("buildBGijCPU [idx=%d]: bgij(%d + %d, %d + %d) [%.15f + %.15fi] != gijTest(%d, %d) [%.15f + %.15fi]\n", idx,
//                      iOffset, i, jOffset, j, bgij(iOffset + i, jOffset + j).real(), bgij(iOffset + i, jOffset + j).imag(),
//                      i, j, gijTest(i, j).real(), gijTest(i, j).imag());
//               exitCompare = true;
//             }
//             if (!comp(bgij(iOffset + kkri + i, jOffset + kkrj + j).real(), gijTest(i, j).real()) ||
//                 !comp(bgij(iOffset + kkri + i, jOffset + kkrj + j).imag(), gijTest(i, j).imag()))
//             {
//               printf("buildBGijCPU : bgij(%d + %d, %d + %d) [%.15f + %.15fi] != gijTest(%d, %d) [%.15f + %.15fi]\n",
//                      iOffset, i + kkri, jOffset, j + kkrj, bgij(iOffset + kkri + i, jOffset + kkrj + j).real(), bgij(iOffset + kkri + i, jOffset + kkrj + j).imag(),
//                      i, j, gijTest(i, j).real(), gijTest(i, j).imag());
//               exitCompare = true;
//             }
//             if (bgij(iOffset + kkri + i, jOffset + j) != 0.0)
//             {
//               printf("buildBGijCPU : bgij(%d + %d, %d + %d) [%.15f + %.15fi] != 0.0\n",
//                      iOffset, i + kkri, jOffset, j, bgij(iOffset + kkri + i, jOffset + j).real(), bgij(iOffset + kkri + i, jOffset + j).imag());
//               exitCompare = true;
//             }
//             if (bgij(iOffset + i, jOffset + kkrj + j) != 0.0)
//             {
//               printf("buildBGijCPU : bgij(%d + %d, %d + %d) [%.15f + %.15fi] != 0.0\n",
//                      iOffset, i, jOffset, j + kkrj, bgij(iOffset + i, jOffset + kkrj + j).real(), bgij(iOffset + i, jOffset + kkrj + j).imag());
//               exitCompare = true;
//             }
//             idx++;
//           }
//       }
//     }
//     if (exitCompare) exit(1);
//   }

// #endif

  smSize = kkrsz_ns * kkrsz_ns * sizeof(deviceFloatComplex);
  threads = 256;

  if (lsms.n_spin_pola ==
      lsms.n_spin_cant)  // non polarized or spin canted
  {
    buildKKRMatrixMultiplyKernelHip_SP<<<blocks, threads>>>(
      devAtom.LIZlmax, devAtom.LIZStoreIdx, devOffsets, kkrsz_ns, ispin,
      lsms.n_spin_pola, lsms.n_spin_cant, iie, d.getBlkSizeTmatStore(),
      d.getTmatStoreLDim(), (deviceFloatComplex *)d.getDevTmatStore_SP(),
      nrmat_ns, (deviceFloatComplex *)devBgij, (deviceFloatComplex *)devM);
  } else {
    buildKKRMatrixMultiplyKernelHip_collinear_SP<<<blocks, threads>>>(
      devAtom.LIZlmax, devAtom.LIZStoreIdx, devOffsets, kkrsz_ns, ispin,
      iie, d.getBlkSizeTmatStore(),
      d.getTmatStoreLDim(), (deviceFloatComplex *)d.getDevTmatStore_SP(),
      nrmat_ns, (deviceFloatComplex *)devBgij, (deviceFloatComplex *)devM);
  }

// #ifdef COMPARE_ORIGINAL
//   Matrix<Complex> mCPU(nrmat_ns, nrmat_ns);
//   Matrix<ComplexF> mGPU(nrmat_ns, nrmat_ns);

//   ret = deviceMemcpy(&mGPU(0, 0), devM, nrmat_ns * nrmat_ns * sizeof(ComplexF),
//              deviceMemcpyDeviceToHost);
//   buildKKRMatrixCPU(lsms, local, atom, ispin, iie, energy, prel, mCPU);

//   for (int i = 0; i < nrmat_ns; i++)
//     for (int j = 0; j < nrmat_ns; j++)
//     {
//       if (!comp(mCPU(i, j).real(), mGPU(i, j).real()) ||
//           !comp(mCPU(i, j).imag(), mGPU(i, j).imag()))
//       {
//         printf("buildBGijCPU : mCPU(%d, %d) [%.15f + %.15fi] != mGPU(%d, %d) [%.15f + %.15fi]\n",
//                i, j, mCPU(i, j).real(), mCPU(i, j).imag(),
//                i, j, mGPU(i, j).real(), mGPU(i, j).imag());
//         exitCompare = true;
//       }
//     }

//   if (exitCompare) exit(1);
// #endif
}

void buildKKRMatrixLMaxDifferentHip_SP(LSMSSystemParameters &lsms,
                                    LocalTypeInfo &local, AtomData &atom,
                                    DeviceStorage &d, DeviceAtom &devAtom,
                                    int ispin, int iie, Complex energy,
                                    Complex prel, ComplexF *devM) {
  int nrmat_ns = lsms.n_spin_cant * atom.nrmat;  // total size of the kkr matrix
  int kkrsz_ns = lsms.n_spin_cant * atom.kkrsz;  // size of t00 block
  bool fullRelativity = false;
  if (lsms.relativity == full) fullRelativity = true;

  ComplexF *devBgij = d.getDevBGijF();

  deviceFloatComplex cuEnergy =
      make_hipFloatComplex(energy.real(), energy.imag());
  deviceFloatComplex cuPrel = make_hipFloatComplex(prel.real(), prel.imag());

  unitMatrixHip<ComplexF>(devM, nrmat_ns, nrmat_ns);
  zeroMatrixHip(devBgij, nrmat_ns, nrmat_ns);

  // calculate Bgij
  // reuse ipvt for offsets
  int *devOffsets = d.getDevIpvt();

  std::vector<int> offsets(devAtom.numLIZ);
  offsets[0] = 0;
  for (int ir = 1; ir < atom.numLIZ; ir++)
    offsets[ir] = offsets[ir - 1] + lsms.n_spin_cant *
                                    (atom.LIZlmax[ir - 1] + 1) *
                                    (atom.LIZlmax[ir - 1] + 1);

  deviceError_t ret = deviceMemcpy(devOffsets, &offsets[0], atom.numLIZ * sizeof(int),
               deviceMemcpyHostToDevice);

  size_t hfnOffset, sinmpOffset, cosmpOffset, plmOffset, dlmOffset;
  size_t smSize = sharedMemoryBGijHip_SP(lsms, &hfnOffset, &sinmpOffset,
                                      &cosmpOffset, &plmOffset, &dlmOffset);
// #ifdef COMPARE_ORIGINAL
//   char *devTestSM;
//   ret = deviceMalloc((void **)&devTestSM, smSize);
// #endif
  int threads = 256;
  dim3 blocks = dim3(devAtom.numLIZ, devAtom.numLIZ, 1);
  buildGijHipKernel_SP<<<blocks, threads, smSize>>>(
      devAtom.LIZPos, devAtom.LIZlmax, DeviceConstants::lofk,
      DeviceConstants::mofk, DeviceConstants::ilp1, DeviceConstants::illp,
      DeviceConstants::cgnt, DeviceConstants::ndlj_illp,
      DeviceConstants::lmaxp1_cgnt, DeviceConstants::ndlj_cgnt, hfnOffset,
      sinmpOffset, cosmpOffset, plmOffset, dlmOffset, cuEnergy, cuPrel,
// #if !defined(COMPARE_ORIGINAL)
      devOffsets, nrmat_ns, (deviceFloatComplex *)devBgij);
// #else
//       devOffsets, nrmat_ns, (deviceFloatComplex *)devBgij, devTestSM);
// #endif
  setBGijHip_SP<<<blocks, threads>>>(fullRelativity, lsms.n_spin_cant,
                                  devAtom.LIZlmax, devOffsets, nrmat_ns,
                                  (deviceFloatComplex *)devBgij);

  smSize = kkrsz_ns * kkrsz_ns * sizeof(deviceFloatComplex);
  threads = 256;

  if (lsms.n_spin_pola ==
            lsms.n_spin_cant)  // non polarized or spin canted
  {
    buildKKRMatrixMultiplyKernelHip_SP<<<blocks, threads>>>(
      devAtom.LIZlmax, devAtom.LIZStoreIdx, devOffsets, kkrsz_ns, ispin,
      lsms.n_spin_pola, lsms.n_spin_cant, iie, d.getBlkSizeTmatStore(),
      d.getTmatStoreLDim(), (deviceFloatComplex *)d.getDevTmatStore_SP(),
      nrmat_ns, (deviceFloatComplex *)devBgij, (deviceFloatComplex *)devM);
  } else {
    buildKKRMatrixMultiplyKernelHip_collinear_SP<<<blocks, threads>>>(
      devAtom.LIZlmax, devAtom.LIZStoreIdx, devOffsets, kkrsz_ns, ispin,
      iie, d.getBlkSizeTmatStore(),
      d.getTmatStoreLDim(), (deviceFloatComplex *)d.getDevTmatStore_SP(),
      nrmat_ns, (deviceFloatComplex *)devBgij, (deviceFloatComplex *)devM);
  }
}

void buildKKRMatrixHip_SP(LSMSSystemParameters &lsms, LocalTypeInfo &local,
                       AtomData &atom, DeviceStorage &devStorage,
                       DeviceAtom &devAtom, int ispin, int iie, Complex energy,
                       Complex prel, ComplexF *devM) {
  // decide between identical lmax and different lmax:

  bool lmaxIdentical = true;

  if (atom.LIZlmax[0] != lsms.maxlmax) {
    lmaxIdentical = false;
    printf("atom.LIZlmax[0] (=%d) != lsms.maxlmax (=%d)\n", atom.LIZlmax[0],
           lsms.maxlmax);
  }
  for (int ir = 0; ir < atom.numLIZ; ir++) {
    if (atom.LIZlmax[ir] != atom.LIZlmax[0]) lmaxIdentical = false;
  }

  if (lmaxIdentical) {
    // printf("lmax identical in buildKKRMatrix\n");
    buildKKRMatrixLMaxIdenticalHip_SP(lsms, local, atom, devStorage, devAtom,
                                   ispin, iie, energy, prel, devM);
  } else {
    // printf("lmax not identical in buildKKRMatrix\n");
    buildKKRMatrixLMaxDifferentHip_SP(lsms, local, atom, devStorage, devAtom,
                                   ispin, iie, energy, prel, devM);
  }
}
