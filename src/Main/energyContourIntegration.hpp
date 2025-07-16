#ifndef LSMS_ENERGYCONTOURINTEGRATION_H
#define LSMS_ENERGYCONTOURINTEGRATION_H

#include "Communication/LSMSCommunication.hpp"
#include "Complex.hpp"
#include "SystemParameters.hpp"

// typedef enum {EnergyGridBox=1, EnergyGridGauss=2} EnergyGridType;

void energyContourIntegration(LSMSCommunication &comm,
                              LSMSSystemParameters &lsms, LocalTypeInfo &local);

extern "C" {
void congauss_(double *ebot, double *etop, double *eibot, DComplex *egrd,
               DComplex *dele1, int *npts, int *nume, double *pi, int *ipepts,
               int *iprint, char *istop, int istop_len);

void congauss_sp_(float *ebot, float *etop, float *eibot, FComplex *egrd,
                  FComplex *dele1, int *npts, int *nume, float *pi, int *ipepts,
                  int *iprint, char *istop, int istop_len);
}

#endif
