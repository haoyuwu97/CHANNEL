#pragma once

#include <channel/config.hpp>
#include <channel/grid.hpp>
#include <channel/kernels.hpp>
#include <channel/stationary.hpp>

namespace channel {

struct VerificationReport {
  bool maxwell_ok = false;
  double maxwell_mismatch = 0.0; // C/m^2

  bool energy_ok = false;
  double dOmega_dVG = 0.0;       // C/m^2 (since Ω per area)
  double Q_gate = 0.0;           // C/m^2
  double energy_rel_error = 0.0;

  bool all_ok() const { return maxwell_ok && energy_ok; }
};

VerificationReport verify_stationary(const ChannelConfig& cfg,
                                     const Grid1D& grid,
                                     const KernelLibrary& kernels,
                                     const StationarySolution& sol);

} // namespace channel
