//
// Created by F.Moitzi on 20.09.2022.
//

#include "SphericalHarmonics.hpp"

int lsms::pvt(int l, int m) { return ((m) + ((l) * ((l) + 1)) / 2); }

int lsms::yvr(int l, int m) { return ((m) + (l) + ((l) * (l))); }

lsms::SphericalHarmonics::SphericalHarmonics(int lmax) : _lmax{lmax} {
  auto sizeP = (lmax + 1) * (lmax + 2) / 2;

  P.resize(sizeP);
  A.resize(sizeP);
  B.resize(sizeP);

  // Precompute coefficients a_l^m and b_l^m for all l <= L, m <= l

  for (int l = 2; l <= lmax; l++) {
    Real ls = l * l;

    Real lm1s = (l - 1) * (l - 1);

    for (int m = 0; m < l - 1; m++) {
      Real ms = m * m;

      A[pvt(l, m)] = std::sqrt((4.0 * ls - 1.0) / (ls - ms));
      B[pvt(l, m)] = -std::sqrt((lm1s - ms) / (4.0 * lm1s - 1.0));
    }
  }
}

void lsms::SphericalHarmonics::computeP(int lmax, Real cos_theta) {
  if (_lmax < lmax) {
    throw std::runtime_error("The arrays are too small!!!");
  }

  const Real sintheta = sqrt(1.0 - cos_theta * cos_theta);
  Real temp = 0.39894228040143267794;  // = sqrt(0.5/M_PI)
  P[pvt(0, 0)] = temp;

  if (lmax > 0) {
    constexpr Real SQRT3 = 1.7320508075688772935;
    P[pvt(1, 0)] = cos_theta * SQRT3 * temp;
    constexpr Real SQRT3DIV2 = -1.2247448713915890491;
    temp = SQRT3DIV2 * sintheta * temp;
    P[pvt(1, 1)] = temp;

    for (int l = 2; l <= lmax; l++) {
      for (int m = 0; m < l - 1; m++) {
        P[pvt(l, m)] = A[pvt(l, m)] * (cos_theta * P[pvt(l - 1, m)] +
                                       B[pvt(l, m)] * P[pvt(l - 2, m)]);
      }
      P[pvt(l, l - 1)] = cos_theta * sqrt(2 * (l - 1) + 3) * temp;
      temp = -sqrt(1.0 + 0.5 / l) * sintheta * temp;
      P[pvt(l, l)] = temp;
    }
  }
}

void lsms::SphericalHarmonics::computeYlm(
    int lmax, std::vector<Real> vec, std::vector<std::complex<Real>> &Ylm) {
  int kmax = (lmax + 1) * (lmax + 1);

  Real q2 = vec[0] * vec[0] + vec[1] * vec[1];
  Real r = std::sqrt(q2 + vec[2] * vec[2]);
  Real q = std::sqrt(q2);

  if (r < TOL) {
    Ylm[0] = 0.28209479177387814347403972578038629292;
    std::fill(Ylm.begin() + 1, Ylm.end(), std::complex<Real>(0.0, 0.0));

  } else {
    using namespace std::complex_literals;

    Real cos_theta = vec[2] / r;

    std::complex<Real> iphi = std::complex<Real>(0, 1) * toReal(atan2(vec[1], vec[0]));

    std::complex<Real> mphi;

    computeP(lmax, cos_theta);

    for (int l = 0; l <= lmax; l++) {
      Ylm[yvr(l, 0)] = P[pvt(l, 0)] * M_SQRT1_2;

      for (int m = 1; m <= l; m++) {
        mphi = iphi * std::complex<Real>(m);

        Ylm[yvr(l, m)] = P[pvt(l, m)] * std::exp(mphi) * toReal(M_SQRT1_2);
        Ylm[yvr(l, -m)] = toReal(std::pow(-1, m)) * std::conj(Ylm[yvr(l, m)]);
      }
    }
  }

  //  for (size_t l = 0; l <= L; l++) {
  //    Y[YVR(l, 0)] = P[PVT(l, 0)] * 0.5 * M_SQRT2;
  //  }
}
