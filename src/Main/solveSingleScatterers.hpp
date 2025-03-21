/* -*- c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */
#include <cmath>
#include <vector>

#include "Complex.hpp"
#include "Misc/Coeficients.hpp"
#include "Misc/Indices.hpp"
#include "PhysicalConstants.hpp"
#include "Real.hpp"
#include "SingleSite/SingleSiteScattering.hpp"

#ifndef LSMS_SOLVE_SINGLE_SCATTERERS
#define LSMS_SOLVE_SINGLE_SCATTERERS

void solveSingleScatterers(
    LSMSSystemParameters &lsms, LocalTypeInfo &local,
    std::vector<Matrix<Real>> &vr, Complex energy,
    std::vector<NonRelativisticSingleScattererSolution> &solution, int iie);

void solveSingleScatterers(
    LSMSSystemParameters &lsms, LocalTypeInfo &local,
    std::vector<Matrix<Real>> &vr, Complex energy,
    std::vector<RelativisticSingleScattererSolution> &solution, int iie);

#endif  // LSMS_SOLVE_SINGLE_SCATTERERS
