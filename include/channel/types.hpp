#pragma once

#include <optional>
#include <string>
#include <vector>

namespace channel {

struct Species {
  std::string name;
  int kappa = 0;                 // valence κ_i
  double q = 0.0;                // charge per particle [C] = κ_i e
  double r_born = 0.0;           // Born radius [m]
  double D = 0.0;                // diffusion coefficient [m^2/s]
  double c_res = 0.0;            // reservoir number density [1/m^3]
  bool mobile = true;
};

struct RedoxParams {
  // W_site = ΔG0 - kBT ln(KX cX) + σP e ψ
  double deltaG0 = 0.0;        // [J]
  double KX = 1.0;             // should satisfy KX*cX dimensionless
  int sigmaP = 1;              // ±1
  std::string counterion = ""; // species label X

  // Kinetics: kon/koff = KX exp[-β(ΔG0 + σP e ψ)]
  // We take koff as a baseline and compute kon locally from the ratio.
  double koff = 0.0;           // [1/s] (0 => equilibrium-only / instantaneous)
};

struct DeviceParams {
  // OECT mapping parameters (optional)
  double Wch = 0.0;  // [m]
  double L = 0.0;    // [m]
  double d = 0.0;    // [m]
  double VD = 0.0;   // [V]

  // mobility model placeholder: µ(ᾱ) = mu0 * exp(mu_alpha * ᾱ)
  double mu0 = 0.0;       // [m^2/V/s]
  double mu_alpha = 0.0;  // dimensionless

  // Map choice: "A" (C*_dyn-based) or "B" (Q-based)
  std::string map = "B";
};

} // namespace channel
