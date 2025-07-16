/* -*- c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */
#ifndef LSMS_SINGLESITESCATTERING_H
#define LSMS_SINGLESITESCATTERING_H

#include "Array3d.hpp"
#include "AtomData.hpp"
#include "Complex.hpp"
#include "Main/SystemParameters.hpp"
#include "Matrix.hpp"

extern "C"
{
void trltog_(int *, int *, Complex *, Complex *, Complex *, Complex *, Complex *);
void gjinv_(Complex *a, int *n, int *nmax, Complex *detl);
void tripmt_(Complex *u, Complex *b, Complex *ust, int *ndi1, int *ind2, int *ndim);
}

class SingleScattererSolution {
 public:
  Complex energy;
  int kkrsz;

  AtomData *atom;
  Matrix<Complex> tmat_g;

  virtual void init(LSMSSystemParameters &lsms, AtomData &a,
                    Complex *tmat_g_store = NULL) {
    printf("Call to 'SingleScattererSolution' base class not allowed!\n");
    exit(1);
  }

  virtual bool isRelativistic() { return false; }
  virtual bool isNonRelativistic() { return false; }
  virtual bool isSpherical() { return false; }
  virtual bool isFullPotential() { return false; }
};

class NonRelativisticSingleScattererSolution : public SingleScattererSolution {
 public:
  NonRelativisticSingleScattererSolution() {}

  NonRelativisticSingleScattererSolution(LSMSSystemParameters &lsms,
                                         AtomData &a,
                                         Complex *tmat_g_stor = NULL) {
    init(lsms, a, tmat_g_stor);
  }

  bool isNonRelativistic() { return true; }
  bool isSpherical() { return true; }

  void init(LSMSSystemParameters &lsms, AtomData &a,
            Complex *tmat_g_store = NULL) {
    atom = &a;
    kkrsz = a.kkrsz;
    matom.resize(a.lmax + 1, 2);
    tmat_l.resize(a.kkrsz, a.kkrsz, 2);
    zlr.resize(a.r_mesh.size(), a.lmax + 1, 2);
    jlr.resize(a.r_mesh.size(), a.lmax + 1, 2);
    if (tmat_g_store != NULL)
      tmat_g.retarget(a.kkrsz * lsms.n_spin_cant, a.kkrsz * lsms.n_spin_pola,
                      tmat_g_store);
    else
      tmat_g.resize(a.kkrsz * lsms.n_spin_cant, a.kkrsz * lsms.n_spin_pola);
  }

  // non relativistic wave functions
  // zlr(ir, l, spin)
  Array3d<Complex> zlr, jlr;
  Matrix<Complex> matom;
  Array3d<Complex> tmat_l;

  Complex ubr[4], ubrd[4];
};

////////////////
// Relativistic
////////////////

// lambdaIndexFromKappaMu
// kappa != 0
// twoMu = 2*mu, i.e. odd integers ..., -7, -5, -3, -1, +1, +3, +5, +7, ...
inline int lambdaIndexFromKappaMu(int kappa, int twoMu) {
  return 2 * kappa * kappa + kappa + (twoMu - 1) / 2;
}

// clebschGordonCoefficientHalf
// Clebsch-Gordon coefficients for j_2 = 1/2

// ms: spin m quantum number. ms=+1/2 -> +1; ms=-1/2 -> -1
inline double clebschGordonCoefficientJ2Half(int kappa, int m, int ms) {
  int l;
  double c, twolp1;
  if (kappa < 0) {
    l = -kappa - 1;
    twolp1 = 2.0 * l + 1.0;
    if (ms < 0) {
      c = std::sqrt((l - m + 1.0) / twolp1);
    } else {
      c = std::sqrt((l + m + 1.0) / twolp1);
    }
  } else {
    l = kappa;
    twolp1 = 2.0 * l + 1.0;
    if (ms < 0) {
      c = std::sqrt((double)(l + m) / twolp1);
    } else {
      c = -std::sqrt((double)(l - m) / twolp1);
    }
  }

  return c;
}

class RelativisticSingleScattererSolution : public SingleScattererSolution {
 public:
  // relativistic wave functions
  static const int nuzp = 2;
  Array3d<Complex> gz, fz, gj, fj;

  std::vector<int> nuz;
  Matrix<int> indz;

  Matrix<Complex> matom;

  RelativisticSingleScattererSolution() {}

  RelativisticSingleScattererSolution(LSMSSystemParameters &lsms, AtomData &a,
                                      Complex *tmat_g_stor = NULL) {
    init(lsms, a, tmat_g_stor);
  }

  bool isRelativistic() { return true; }
  bool isSpherical() { return true; }

  void init(LSMSSystemParameters &lsms, AtomData &a,
            Complex *tmat_g_store = NULL) {
    atom = &a;
    kkrsz = a.kkrsz;
    matom.resize(a.lmax + 1, 2);

    if (tmat_g_store != NULL)
      tmat_g.retarget(a.kkrsz * lsms.n_spin_cant, a.kkrsz * lsms.n_spin_pola,
                      tmat_g_store);
    else
      tmat_g.resize(a.kkrsz * lsms.n_spin_cant, a.kkrsz * lsms.n_spin_pola);

    gz.resize(a.r_mesh.size(), nuzp, 2 * a.kkrsz);
    fz.resize(a.r_mesh.size(), nuzp, 2 * a.kkrsz);
    gj.resize(a.r_mesh.size(), nuzp, 2 * a.kkrsz);
    fj.resize(a.r_mesh.size(), nuzp, 2 * a.kkrsz);

    nuz.resize(2 * a.kkrsz);
    indz.resize(nuzp, 2 * a.kkrsz);

    // printf("relativistiv wavefunction not implemented yet!\n");
    // exit(1);
  }
};

// Wave functions are expanded in free space spherical solutions using the phase
// method as phi_L(r) = sum_L' sineMatrix(r,L',L) * N_L'(r) -
// cosineMatrix(r,L',L) * J_L'(r)

class NewRelativisticSingleScattererSolution : public SingleScattererSolution {
 public:
  Array3d<Complex> sineMatrixRegular,
      cosineMatrixRegular;  // sineMatrix(ir, Lambda', Lambda)
  Array3d<Complex> sineMatrixIrregular,
      cosineMatrixIrregular;  // sineMatrix(ir, Lambda', Lambda)

  bool isRelativistic() { return true; }
  bool isSpherical() { return true; }
};

class FullPotentialRelativisticSingleScattererSolution
    : public SingleScattererSolution {
 public:
  Array3d<Complex> sineMatrixRegular,
      cosineMatrixRegular;  // sineMatrix(ir, Lambda', Lambda)
  Array3d<Complex> sineMatrixIrregular,
      cosineMatrixIrregular;  // sineMatrix(ir, Lambda', Lambda)

  bool isRelativistic() { return true; }
  bool isFullPotential() { return true; }
};

extern "C" {
void single_site_tmat_(int *nrel_rel, int *n_spin_cant, int *is,
                       int *n_spin_pola, int *mtasa, double *rws, int *nrelv,
                       double *clight, int *lmax, int *kkrsz, DComplex *energy,
                       DComplex *prel, DComplex *pnrel, double *vr, double *h,
                       int *jmt, int *jws, double *r_mesh, DComplex *tmat_l,
                       DComplex *tmat_g, DComplex *matom, DComplex *zlr,
                       DComplex *jlr, DComplex *gz, DComplex *fz, DComplex *gj,
                       DComplex *fj, int *nuz, int *indz, DComplex *ubr,
                       DComplex *ubrd, DComplex *dmat, DComplex *dmatp,
                       double *r_sph, int *iprint, const char *istop);

void single_site_tmat_sp_(int *nrel_rel, int *n_spin_cant, int *is,
                       int *n_spin_pola, int *mtasa, float *rws, int *nrelv,
                       float *clight, int *lmax, int *kkrsz, FComplex *energy,
                       FComplex *prel, FComplex *pnrel, float *vr, float *h,
                       int *jmt, int *jws, float *r_mesh, FComplex *tmat_l,
                       FComplex *tmat_g, FComplex *matom, FComplex *zlr,
                       FComplex *jlr, FComplex *gz, FComplex *fz, FComplex *gj,
                       FComplex *fj, int *nuz, int *indz, FComplex *ubr,
                       FComplex *ubrd, FComplex *dmat, FComplex *dmatp,
                       float *r_sph, int *iprint, const char *istop);

void single_scatterer_nonrel_(int *nrelv, double *clight, int *lmax, int *kkrsz,
                              DComplex *energy, DComplex *prel, DComplex *pnrel,
                              double *vr, double *r_mesh, double *h, int *jmt,
                              int *jws, DComplex *tmat_l, DComplex *matom,
                              DComplex *zlr, DComplex *jlr, double *r_sph,
                              int *iprpts, int *iprint, char *istop,
                              int istop_len);

void single_scatterer_nonrel_sp_(int *nrelv, float *clight, int *lmax, int *kkrsz,
                              FComplex *energy, FComplex *prel, FComplex *pnrel,
                              float *vr, float *r_mesh, float *h, int *jmt,
                              int *jws, FComplex *tmat_l, FComplex *matom,
                              FComplex *zlr, FComplex *jlr, float *r_sph,
                              int *iprpts, int *iprint, char *istop,
                              int istop_len);

void single_scatterer_rel_(
    DComplex *ce, DComplex *psq, int *lmax, int *kmymax,
    int *idpot,  // idpot=identifies vauum empty sphere: idpot=0 -> Vacuum else
                 // an atomic site (can be set to Z!)
    double *v0, double *vr, double *br, double *bopr, double *dx, int *ns,
    double *rs, DComplex *tminv, DComplex *gz, DComplex *fz, DComplex *gj,
    DComplex *fj, int *nuz, int *indz, int *iflag, double *socsc, int *iprpts,
    int *iprint, const char *istop, int istop_len);

void single_scatterer_rel_sp_(
    FComplex *ce, FComplex *psq, int *lmax, int *kmymax,
    int *idpot,  // idpot=identifies vauum empty sphere: idpot=0 -> Vacuum else
                 // an atomic site (can be set to Z!)
    float *v0, float *vr, float *br, float *bopr, float *dx, int *ns,
    float *rs, FComplex *tminv, FComplex *gz, FComplex *fz, FComplex *gj,
    FComplex *fj, int *nuz, int *indz, int *iflag, float *socsc, int *iprpts,
    int *iprint, const char *istop, int istop_len);    
}

void calculateSingleScattererSolution(
    LSMSSystemParameters &lsms, AtomData &atom, Matrix<Real> &vr,
    Complex energy, Complex prel, Complex pnrel,
    NonRelativisticSingleScattererSolution &solution);

void calculateScatteringSolutions(
    LSMSSystemParameters &lsms, std::vector<AtomData> &atom, Complex energy,
    Complex prel, Complex pnrel,
    std::vector<NonRelativisticSingleScattererSolution> &solution);

void calculateSingleScattererSolution(
    LSMSSystemParameters &lsms, AtomData &atom, Matrix<Real> &vr,
    Complex energy, RelativisticSingleScattererSolution &solution);

#endif
