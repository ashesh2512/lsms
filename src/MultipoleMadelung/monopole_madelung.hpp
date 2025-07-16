//
// Created by F.Moitzi on 06.01.2022.
//

#ifndef LSMS_MONO_MADELUNG_HPP
#define LSMS_MONO_MADELUNG_HPP

#include <complex>
#include <vector>

#include "Real.hpp"
#include "common.hpp"

namespace lsms {

constexpr auto EPSI = 1e-14;

/**
 * Calculate the scaling factor to get a balanced number
 * of reciprocal and real-space lattice vectors
 */
Real scaling_factor(const matrix<Real> &bravais, int lmax,
                      int max_iter = 1000000, Real fstep = 0.2);

/**
 * Number of lattice vectors
 */
int num_latt_vectors(const matrix<Real> &brav, Real cut,
                     const std::vector<int> &nm);

/**
 * Get radius of truncation sphere
 */
Real rs_trunc_radius(const matrix<Real> &brav, int lmax, Real eta,
                       const std::vector<int> &nm);

Real kn_trunc_radius(const matrix<Real> &brav, int lmax, Real eta,
                       const std::vector<int> &nm);

/**
 * Get size of lattice multiplications
 */
std::vector<int> real_space_multiplication(const matrix<Real> &brav, int lmax,
                                           Real eta);

/**
 * Get size of reciprocal lattice multiplications
 */
std::vector<int> reciprocal_space_multiplication(const matrix<Real> &brav,
                                                 int lmax, Real eta);

/**
 * Calculate the `\eta` factor
 */
Real calculate_eta(matrix<Real> &brav);

}  // namespace lsms

#endif  // LSMS_MONO_MADELUNG_HPP
