//
// Created by F.Moitzi on 10.01.2023.
//

#ifndef LSMS_SRC_CORE_ATOMIC_DFT_HPP_
#define LSMS_SRC_CORE_ATOMIC_DFT_HPP_

#include <vector>
#include <memory>

#include "SingleSite/AtomData.hpp"
#include "XCLDA.hpp"

namespace lsms {

constexpr int BROYDEN_MAX_STEP = 6;
constexpr int BROYDEN_ITER_RESET = 25;
constexpr Real BROYDEN_W0 = 0.01;

class AtomicDFT {

 private:

  int max_iter;
  int max_eig_iter;
  Real e_eig_tol;
  Real e_tol;
  int iprint;
  Real alpha;

  std::unique_ptr<lsms::XCBase> xc;

 public:

  explicit AtomicDFT(
      std::vector<int> functional = {0, 0, 0},
      int max_iter = 100,
      int max_eig_iter = 100,
      Real e_tol = 1.0e-10,
      Real e_eig_tol = 1.0e-12, int iprint = 0, Real alpha = 0.05
  );

  std::tuple<std::vector<Real>, Real> solve(int Z,
               const std::vector<Real> &r_mesh,
               Real h,
               int N,
               Matrix<Real> &density,
               std::vector<Real> &tot_potential
  );

};

void generate_starting_potential(std::vector<Real> &potential,
                                 const std::vector<Real> &r_mesh,
                                 int core_charge,
                                 std::size_t end
);

Real generate_density(
    const std::vector<Real> &r_mesh,
    Real h,
    std::size_t end,
    const std::vector<Real> &potential,
    int core_charge,
    const std::vector<int> &n,
    const std::vector<int> &l,
    const std::vector<int> &spin,
    const std::vector<int> &kappa,
    const std::vector<Real> &occupation,
    std::vector<Real> &e_eig,
    Matrix<Real> &density,
    int max_eig_iter,
    Real eig_tol
);

Real total_energy(const std::vector<Real> &r_mesh,
                    Matrix<Real> &rho,
                    std::vector<Real> &v_hartree,
                    Matrix<Real> &e_xc,
                    std::vector<Real> &v_pot,
                    Real E_band,
                    Real Z,
                    int N
);

} // lsms

#endif //LSMS_SRC_CORE_ATOMIC_DFT_HPP_
