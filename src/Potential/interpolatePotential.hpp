#ifndef INTERPOLATE_POTENTIAL_HPP
#define INTERPOLATE_POTENTIAL_HPP

#include "Main/SystemParameters.hpp"
#include "SingleSite/AtomData.hpp"

extern "C" {
void interp_(double *r, double *f, int *nr, double *rs, double *ps, double *dps,
             int *deriv);
void interp_sp_(float *r, float *f, int *nr, float *rs, float *ps, float *dps,
                int *deriv);
}

void interpolatePotential(LSMSSystemParameters &lsms, AtomData &atom);

#endif
