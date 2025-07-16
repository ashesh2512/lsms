#ifndef LSMS_PHYSICALCONSTANTS_H
#define LSMS_PHYSICALCONSTANTS_H

#include "Real.hpp"

// fine-structure constant alpha
constexpr Real alphaInverse = 137.035999206;
constexpr Real alpha = 1/alphaInverse;
constexpr Real cphot = 2 * alphaInverse; // 274.072;
constexpr Real c2inv = 1/(cphot*cphot);

// Joule in Rydberg
constexpr Real convertJouleToRydberg = 1.380649e-23;
// Rydberg in eV
constexpr Real convertRydbergToeV = 13.605698066;

// Boltzman constant in Ry/K
// defined by SI to be 1.380649e-23 J/K

// const double kBoltzmann = 6.3336e-6;
constexpr Real kBoltzmann = 1.380649e-23 * 4.5874208973812E+17;

constexpr Real convertKtoRydberg = 6.33361706838587e-06;

constexpr long double PI = 3.141592653589793238462643383279502884L;

#endif
