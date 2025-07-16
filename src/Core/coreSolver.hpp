/* -*- c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */
#ifndef LSMS_CORESOLVER_HPP
#define LSMS_CORESOLVER_HPP

#include "Communication/LSMSCommunication.hpp"
#include "Main/SystemParameters.hpp"

extern "C" {
void deepst_(int *nqn, int *lqn, int *kqn, double *en, double *rv, double *r,
             double *rf, double *h, double *z, double *c, int *nitmax, double *tol,
             int *nws, int *nlast, int *iter, int *iprpts, int *ipdeq);
void deepst_sp_(int *nqn, int *lqn, int *kqn, float *en, float *rv, float *r,
                float *rf, float *h, float *z, float *c, int *nitmax, float *tol,
                int *nws, int *nlast, int *iter, int *iprpts, int *ipdeq);

void semcst_(int *nqn, int *lqn, int *kqn, double *en, double *rv, double *r,
             double *rf, double *h, double *z, double *c, int *nitmax, double *tol,
             int *nmt, int *nws, int *nlast, int *iter, int *iprpts,
             int *ipdeq);
void semcst_sp_(int *nqn, int *lqn, int *kqn, float *en, float *rv, float *r,
                float *rf, float *h, float *z, float *c, int *nitmax, float *tol,
                int *nmt, int *nws, int *nlast, int *iter, int *iprpts,
                int *ipdeq);
             
void newint_(int *nr, double *r, double *f, double *g, int *ip0);
void newint_sp_(int *nr, float *r, float *f, float *g, int *ip0);
}

void getCoreStates(LSMSSystemParameters &lsms, AtomData &atom);

namespace lsms {

const static int NITMAX = 50;
const static int IPDEQ = 5;

#if SP
const static Real TOL = 1.0e-05; // TOLERANCE FOR SP
#else
const static Real TOL = 1.0e-10;
#endif

void getNewCoreStates(LSMSSystemParameters &lsms, AtomData &atom);

}  // namespace lsms

#endif
