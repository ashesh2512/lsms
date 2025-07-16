#include "bulirschStoerIntegrator.hpp"

#include <cmath>
#include <complex>

/*
      subroutine bulirsch_stoer(y,dy,nv,x0,dx,
     &                          rn,ecomp,vr,eunit,b,allp1,nextrp)
*/

extern "C" {
void dfv_m_(double *r, double *y, double *dy, int *nv, int *n, double *rn,
            std::complex<double> *ecomp, double *vr, double *eunit, double *b,
            double *allp1);

void dfv_m_sp_(float *r, float *y, float *dy, int *nv, int *n, float *rn,
            std::complex<float> *ecomp, float *vr, float *eunit, float *b,
            float *allp1);

void bulirsch_stoer_integrator_(Real *xFrom, Real *xTo, Real *y,
                                Real *dy, int *nv, Real *rn,
                                Complex *ecomp, Real *vr,
                                Real *eunit, Real *b, Real *allp1,
                                int *nextrp);
}

// specialization to use the dfv_m subroutine to calculate the right hand side
// of the radial ODE
inline void modifiedMidpoint_dfv(Real x0, Real x1, Real *y0, Real *y1,
                                 int n, int steps, int *nv, int *nextrp,
                                 Real *rn, Complex *ecomp,
                                 Real *vr, Real *eunit, Real *b,
                                 Real *allp1) {
  Real y2[n], dy[n];
  Real h = (x1 - x0) / Real(steps);
  Real h2 = Real(2) * h;

  Real x = x0;

  // rhs(x,y0,dy);
#if SP
  dfv_m_sp_(&x, y0, dy, nv, nextrp, rn, ecomp, vr, eunit, b, allp1);
#else
  dfv_m_(&x, y0, dy, nv, nextrp, rn, ecomp, vr, eunit, b, allp1);
#endif

  for (int i = 0; i < n; i++) {
    y1[i] = y0[i];
    y2[i] = y0[i] + h * dy[i];
  }
  for (int j = 1; j < steps; j++) {
    x += h;
    // rhs(x,y2,dy);
#if SP
    dfv_m_sp_(&x, y2, dy, nv, nextrp, rn, ecomp, vr, eunit, b, allp1);
#else
    dfv_m_(&x, y2, dy, nv, nextrp, rn, ecomp, vr, eunit, b, allp1);
#endif
    for (int i = 0; i < n; i++) {
      Real temp = y1[i] + h2 * dy[i];
      y1[i] = y2[i];
      y2[i] = temp;
    }
  }
  x += h;
  // rhs(x,y2,dy);
#if SP
  dfv_m_sp_(&x, y2, dy, nv, nextrp, rn, ecomp, vr, eunit, b, allp1);
#else
  dfv_m_(&x, y2, dy, nv, nextrp, rn, ecomp, vr, eunit, b, allp1);
#endif
  for (int i = 0; i < n; i++) {
    y1[i] = Real(0.5) * (y1[i] + y2[i] + h * dy[i]);
  }
}

inline int bulirschStoerIntegrator_dfv(Real x0, Real x1, Real *y0,
                                       Real *y1, int n, int *nv, int *nextrp,
                                       Real *rn, std::complex<Real> *ecomp,
                                       Real *vr, Real *eunit, Real *b,
                                       Real *allp1, double eps = 1.0e-12) {
  const int nSteps = 11;
  const int stepTarget = 7;
  const Real stepGrowth = 1.25;
  const Real stepShrink = 0.8;
  // Bulirsch Stoer sequence (n_j = 2 n_{j-2})
  const int stepSequence[] = {2,  4,  6,   8,   12,  16,  24,  32,  48,
                              64, 96, 128, 192, 256, 384, 512, 768, 1024};
  // Deuflhard sequence (n_j = 2 j)
  // const int stepSequence[]={2,4,6,8,10,12,14,16,18,20,22,24,26,28,30};
  // Romberg sequence (n_j = 2^j)
  // const int
  // stepSequence[]={2,4,8,16,32,64,128,256,512,1024,2048,4096,8192,16384,32768,65536};
  // Harmonic sequence (n_j = j)
  // const int stepSequence[]={1,2,3,4,5,6,7,8,9,10,11,12,13,14}

  Real xEpsilon = eps * std::abs(x0); // TOLERANCE FOR SP
  Real err;
  Real y[n];

  Real x = x0;
  Real h0 = x1 - x0;
  Real h = h0;
  for (int i = 0; i < n; i++) y[i] = y0[i];

  ExtrapolatorT0<Real, Real> extrapolation(n, nSteps + 1);

  do {
    // fill extrapolation Tableau
    int step = 0;
    Real err = HUGE_VAL;
    extrapolation.reset();
    while (step < nSteps && err > eps) {
      modifiedMidpoint_dfv(x, x + h, y, y1, n, stepSequence[step], nv, nextrp,
                           rn, ecomp, vr, eunit, b, allp1);
      extrapolation.add(h / Real(stepSequence[step]), y1);
      err = extrapolation.relativeError();
      step++;
    }
    // printf("step %d err %g\n",step,err);
    if (err < eps) {
      x = x + h;
      if (step < stepTarget)
        h *= stepGrowth;
      else
        h *= Real(stepSequence[stepTarget]) / Real(stepSequence[step]);
      h = std::copysign(std::min(std::abs(h), std::abs(x1 - x)), h0);
      extrapolation.estimate(y);
    } else {
      // the error is to large -> reduce the interval
      h = Real(0.25) * h;
      if (h + x == x) {
        printf("Bulirsch-Stoer didn't converge to %g (error=%g)\n", eps, err);
        extrapolation.estimate(y1);
        return 1;
      }
    }
  } while (std::abs(x - x1) > xEpsilon);

  extrapolation.estimate(y1);
  return 0;
}

void bulirsch_stoer_integrator_(Real *xFrom, Real *xTo, Real *y,
                                Real *dy, int *nv, Real *rn,
                                Complex *ecomp, Real *vr,
                                Real *eunit, Real *b, Real *allp1,
                                int *nextrp) {
  Real y1[*nv];
  //  bulirschStoerIntegrator<Real,Real>(*xFrom, *xTo, y, y1, *nv,
  //                                  [&nv,&rn,&ecomp,&vr,&eunit,&b,&allp1,&nextrp](Real
  //                                  x, Real *y, Real *dy)
  //                                  {
  //                                    dfv_m_(&x, y, dy, nv, nextrp, rn, ecomp,
  //                                    vr, eunit, b, allp1);
  //                                  },
  //                                  1.0e-12);
  bulirschStoerIntegrator_dfv(*xFrom, *xTo, y, y1, *nv, nv, nextrp, rn, ecomp,
                              vr, eunit, b, allp1, 1.0e-12);

  for (int i = 0; i < *nv; i++) y[i] = y1[i];

#if SP
  dfv_m_sp_(xTo, y, dy, nv, nextrp, rn, ecomp, vr, eunit, b, allp1);
#else
  dfv_m_(xTo, y, dy, nv, nextrp, rn, ecomp, vr, eunit, b, allp1);
#endif
}
