/* -*- c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */
#ifndef LSMS_CORESTATES_H
#define LSMS_CORESTATES_H

#include "Main/SystemParameters.hpp"
#include "Communication/LSMSCommunication.hpp"

namespace lsms {

void calculateCoreStates(LSMSCommunication &comm, LSMSSystemParameters &lsms, LocalTypeInfo &local);

void calculateCoreStates(LSMSCommunication &comm, LSMSSystemParameters &lsms, AlloyAtomBank &alloyBank);

/*
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     subroutine getcor(IN n_spin_pola,IN mtasa,
!    >                  IN jmt,IN jws,IN r_mesh,IN h,IN xstart,IN vr,
!    >                  IN numc,IN nc,IN lc,IN kc,INOUT ec,
!    >                  IN ztotss,IN zsemss,IN zcorss,
!    >                  OUT ecorv,OUT esemv,OUT corden,OUT semcor,
!    >                  IN nrelc,
!    >                  OUT qcpsc_mt,OUT qcpsc_ws,OUT mcpsc_mt,OUT mcpsc_ws,
!    >                  IN iprpts, IN ipcore,
!    >                  IN iprint,IN istop)
c     ================================================================
*/

}

extern "C" {
void getcor_(int *n_spin_pola, int *mtasa, int *jmt, int *jws, double *r_mesh,
             double *h, double *xstart, double *vr, int *numc, int *nc, int *lc,
             int *kc, double *ec, double *ztotss, double *zsemss,
             double *zcorss, double *ecorv, double *esemv, double *corden,
             double *semcor, int *nrelc, double *qcpsc_mt, double *qcpsc_ws,
             double *mcpsc_mt, double *mcpsc_ws, int *iprpts, int *ipcore,
             int *iprint, char *istop, int itop_length);

void getcor_sp_(int *n_spin_pola, int *mtasa, int *jmt, int *jws, float *r_mesh,
               float *h, float *xstart, float *vr, int *numc, int *nc, int *lc,
               int *kc, float *ec, float *ztotss, float *zsemss,
               float *zcorss, float *ecorv, float *esemv, float *corden,
               float *semcor, int *nrelc, float *qcpsc_mt, float *qcpsc_ws,
               float *mcpsc_mt, float *mcpsc_ws, int *iprpts, int *ipcore,
               int *iprint, char *istop, int itop_length);
}

#endif
