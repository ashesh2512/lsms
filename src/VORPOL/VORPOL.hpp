/* -*- c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */
#ifndef LSMS_VORPOL_H
#define LSMS_VORPOL_H

#include "Array3d.hpp"
#include "Complex.hpp"
#include "Matrix.hpp"
#include "Real.hpp"

class VoronoiPolyhedra {
 public:
  /*
  // VORPOL data [this is only used from routines imported from LSMS_1.9 ..
  VORPOL and replaces the common block /common_step/ integer    ncrit real*8
  gwwylm(ipngaussr*(iprcrit-1)) real*8     grwylm(ipngaussr*(iprcrit-1)) include
  'atom_param.h' complex*16 wylm((2*iplmax+1)*(iplmax+1),ipngaussr,iprcrit-1)

  c This common block is defined in subroutine setup_vorpol. It
  c contains the step function and the r-mesh for interstitial integration
  c ncrit = number of critical r points
  c grwylm = array containing the r-mesh for interstitial integration
  c gwwylm = array containing the weights for interstitial integration
  c wylm = array containing the step function for each r-mesh point
   common /common_step/ wylm,gwwylm,grwylm,ncrit
  */
  int ncrit;
  double rInscribedSphere;  // inscribed sphere
  double omegaInt;          // interstitial volume
  DComplex dipint[6];
  Array3d<Complex> wylm;
  Matrix<Real> gwwylm, grwylm;
};

extern "C" {
void setup_vorpol_(int *my_atom, int *num_atoms, double *atom_position_1,
                   double *atom_position_2, double *atom_position_3,
                   double *system_bravais, int *lmax, double *clm, int *ngaussq,
                   int *ngaussr, double *rmt, double *omegint, DComplex *dipint,
                   double *rad, int *ipvp, int *ipnode, int *ipcorn, int *ipedge,
                   int *iprcrit, double *gwwylm, double *grwylm, int *ncrit,
                   DComplex *wylm, double *rcirc, int *iprint, char *istop,
                   int istopl_len);    

void setup_vorpol_vplane_(double *vp, int *nvplane, int *lmax, double *clm,
                          int *ngaussq, int *ngaussr, double *rmt, double *omegint,
                          DComplex *dipint, int *ipvp, int *ipnode, int *ipcorn,
                          int *ipedge, int *iprcrit, double *gwwylm, double *grwylm,
                          int *ncrit, DComplex *wylm, double *rcirc, int *iprint,
                          char *istop, int istopl_len);

void setup_boundary_(int *my_atom, int *num_atoms, double *atom_position_1,
                     double *atom_position_2, double *atom_position_3,
                     double *bravais_1, double *bravais_2, double *bravais_3,
                     double *vplanes, int *ipvp, int *nvplanes, double *rad);

void setup_boundary_cluster_(int *i_seed, int *num_seeds, double *atom_position_1,
                             double *atom_position_2, double *atom_position_3,
                             double *vplanes, int *ipvp, int *nvplanes,
                             double *rad);                        
}

#endif
