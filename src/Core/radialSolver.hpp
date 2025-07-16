
#ifndef LSMS_RADIAL_SOLVER_H
#define LSMS_RADIAL_SOLVER_H

#include <cfloat>
#include <cstdlib>
#include <cmath>

#include "Real.hpp"

namespace lsms {

constexpr int ERROR_WRONG_NUMBER_NODES = -3;
constexpr int ERROR_STILL_BISECTION = -2;
constexpr int ERROR_NOT_REACHED_ACCURACY = -1;

constexpr Real U_LIM = 1e5;
constexpr Real L_LIM = 1e-5;
constexpr Real C = 137.0359991; // Speed of light in atomic units (inverse \f$\alpha\f$)
constexpr Real C_SQ = 18778.86504933520081;

template<typename T, typename D>
inline void rel_rad_func(T r, T rp,
                         T E, T vpot,
                         D k_qn, T p,
                         T q, T val[2]) {
  val[0] = (-k_qn / r * p + ((E - vpot) / C + 2 * C) * q) * rp;
  val[1] = (-((E - vpot) / C) * p + k_qn / r * q) * rp;
}

template<typename T, typename D>
inline void nrel_rad_func(T r, T rp,
                          T E, T vpot,
                          D l, T p, T q,
                          T val[2]) {
  val[0] = q * rp;
  val[1] = toReal(2.0) * ((vpot - E + (l * l + l) / (2 * r * r)) * p) * rp;
}

template<typename T, typename D, typename V>
inline V calculate_ctp(T E, D l, const T *R,
                       const T *pot, V end) {
  auto idx = end;

  for (idx = end - 1; idx >= 0; idx--) {
    if ((E - pot[idx] - (l * l + l) / (toReal(2.0) * R[idx] * R[idx])) > 0) {
      return idx;
    }
  }

  return end;
}

template<typename T, typename D, typename V>
inline int count_sign_sign_switch(const T *P, D start,
                                  V end) {
  auto count = 0;

  for (auto idx = start; idx < end - 1; idx++) {
    if (P[idx] * P[idx + 1] < 0) {
      count += 1;
    }
  }

  return count;
}

template<typename T, typename D>
inline D last_valid_point(const T *P, D end) {
  auto i = end - 1;

  for (i = 0; i < end; i++) {
    if (std::fabs(P[i]) >= lsms::U_LIM) {
      break;
    }
  }

  return i;

}

void nrel_out(Real Z, Real energy, int l_qn, Real *p_array, Real *q_array,
              const Real *__restrict__ r_mesh, Real h,
              const Real *__restrict__ pot, std::size_t end, int sign_switch);

void nrel_in(Real energy, int l_qn, Real *p_array, Real *q_array,
             const Real *__restrict__ r_mesh, Real h,
             const Real *__restrict__ pot, std::size_t end, Real &p_last,
             Real &q_last, std::size_t stop, std::size_t &imax);

void rel_in(Real energy, int k_qn, Real *p_array, Real *q_array,
            const Real *__restrict__ r_mesh, Real h,
            const Real *__restrict__ pot, std::size_t end, Real &p_last,
            Real &q_last, std::size_t stop, std::size_t &imax);

void rel_out(Real energy, int k_qn, Real *p_array, Real *q_array,
             const Real *__restrict__ r_mesh, Real h,
             const Real *__restrict__ pot, std::size_t end, int sign_switch);

inline Real rel_integration(const Real *__restrict__ p_array,
                              const Real *__restrict__ q_array,
                              const Real *__restrict__ r_mesh,
                              Real h, std::size_t end);

inline Real nrel_integration(const Real *__restrict__ P,
                               const Real *__restrict__ r_mesh,
                               Real h,
                               std::size_t end);

Real nonrel_eigenenergies(Real Z, int n_qn, int l_qn, int k_qn, Real *dens,
                            Real *p_array, Real *q_array, const Real *__restrict__ r_mesh,
                            Real h,
                            const Real *__restrict__ pot, std::size_t end,
                            Real energy_rel_tol, int max_iter, Real e_init,
                            int &converged, Real &delta_energy);

Real rel_eigenenergies(
    Real Z, int n_qn, int l_qn, int k_qn, Real *dens, Real *p_array, Real *q_array,
    const Real *__restrict__ r_mesh, Real h,
    const Real *__restrict__ pot, std::size_t end, Real energy_rel_tol, int max_iter,
    Real e_init, int &converged, Real &delta_energy);

Real rel_energy_start(int n, int k_qn, Real Z);

Real nrel_energy_start(int n, Real Z);

}  // namespace lsms

#endif  // LSMS_RADIAL_SOLVER_H
