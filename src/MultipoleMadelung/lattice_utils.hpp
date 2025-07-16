//
// Created by F.Moitzi on 21.12.2021.
//

#ifndef LSMS_LATTICE_HPP
#define LSMS_LATTICE_HPP

#include <complex>
#include <tuple>
#include <vector>

#include "Real.hpp"
#include "common.hpp"

namespace lsms {

/**
 * Create all vectors in certain cutoff with a certain repeotition
 */
__attribute__((unused)) matrix<Real> create_lattice(
    const matrix<Real> &brav, Real cutoff, const std::vector<int> &nm,
    int size);

/**
 * Create all vectors in certain cutoff with a certain repeotition
 */
std::tuple<matrix<Real>, std::vector<Real>> create_lattice_and_sq(
    matrix<Real> &brav, Real cutoff, const std::vector<int> &nm, int size);

/**
 *  inserts a vector in a list of vectors such that they are in
 *  order of increasing length.
 */
void insert_ordered(matrix<Real> &latt_vec, std::vector<Real> &latt_vec_sq,
                    int len, std::vector<Real> &vec, Real &v_sq);

/**
 * Lattice volumes
 */
Real omega(matrix<Real> &bravais);

/**
 * Calculate reciprocal lattice
 */
void reciprocal_lattice(matrix<Real> &bravais,
                        matrix<Real> &reciprocal_bravais, Real &scale);

}  // namespace lsms

#endif  // LSMS_LATTICE_HPP
