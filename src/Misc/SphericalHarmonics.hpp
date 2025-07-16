//
// Created by F.Moitzi on 20.09.2022.
//

#ifndef LSMS_SPHERICALHARMONICS_HPP
#define LSMS_SPHERICALHARMONICS_HPP

#include <cmath>
#include <complex>
#include <vector>

#include "Real.hpp"

namespace lsms {

constexpr static Real TOL = 0.5 * 1.0e-12;

inline int pvt(int l, int m);

inline int yvr(int l, int m);

class SphericalHarmonics {
 public:
  explicit SphericalHarmonics(int lmax);

  /* Compute an entire set of Y_{l,m}(\theta,\phi) and store in array Y */
  void computeYlm(int lmax, std::vector<Real> vec,
                  std::vector<std::complex<Real>> &Ylm);

 private:
  int _lmax;

  std::vector<Real> A;
  std::vector<Real> B;
  std::vector<Real> P;

  /**
   * Compute an entire set of P_l^m(x) and store in the array P
   *
   * @param lmax
   * @param X
   */
  void computeP(int lmax, Real X);
};

}  // namespace lsms

#endif  // LSMS_SPHERICALHARMONICS_HPP
