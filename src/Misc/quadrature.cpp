
#include "Real.hpp"
#include "calculateGaussLegendrePoints.hpp"

extern "C" {
void gauss_legendre_points_sp_(float *x1, float *x2, float *x, float *w,
                               int *n) {
  calculateGaussLegendrePoints<float>(x, w, *n, *x1, *x2);
}

void gauss_legendre_points_dp_(double *x1, double *x2, double *x, double *w,
                               int *n) {
  calculateGaussLegendrePoints<double>(x, w, *n, *x1, *x2);
}
}
