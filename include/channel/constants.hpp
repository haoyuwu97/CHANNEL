#pragma once

namespace channel {

// Fundamental constants (CODATA 2018/2019 level; good enough for this workflow).
constexpr double kBoltzmann = 1.380649e-23;      // J/K
constexpr double eCharge    = 1.602176634e-19;   // C
constexpr double epsilon0   = 8.8541878128e-12;  // F/m
constexpr double pi         = 3.141592653589793238462643383279502884;

inline double beta(double T) { return 1.0 / (kBoltzmann * T); }

} // namespace channel
