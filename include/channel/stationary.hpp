#pragma once

#include <channel/config.hpp>
#include <channel/grid.hpp>
#include <channel/kernels.hpp>
#include <channel/observables.hpp>
#include <channel/types.hpp>

#include <optional>
#include <string>
#include <vector>

namespace channel {

struct StationarySolution {
  double VG = 0.0;
  int iterations = 0;
  bool converged = false;

  std::vector<double> psi;                   // [nz]
  std::vector<double> alpha;                 // [nz]
  std::vector<std::vector<double>> c;        // [nspecies][nz]

  KernelState K;

  ChargeSummary charge;
  double Omega = 0.0;

  // Optional capacitance estimate around VG (finite difference).
  bool has_capacitance = false;
  double C_eq = 0.0; // [F/m^2] per-area capacitance
};

StationarySolution solve_stationary(const ChannelConfig& cfg,
                                    const Grid1D& grid,
                                    const KernelLibrary& kernels,
                                    double VG_override,
                                    bool compute_capacitance);

void write_stationary_outputs(const ChannelConfig& cfg,
                              const Grid1D& grid,
                              const std::vector<Species>& species,
                              const StationarySolution& sol,
                              const std::string& subdir = "");

} // namespace channel
