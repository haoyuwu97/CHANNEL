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

struct DynamicSolution {
  std::vector<double> t;         // [s]
  std::vector<double> VG;        // [V]
  std::vector<double> Q_gate;    // [C/m^2]
  std::vector<double> Q_vol;     // [C/m^2]
  std::vector<double> alpha_bar; // [-]
  std::vector<double> Omega;     // [J/m^2] (for diagnostics)

  // Optional mapped device current
  std::vector<double> ID;        // [A] if device enabled (else empty)

  // Final fields
  std::vector<double> psi;
  std::vector<double> alpha;
  std::vector<std::vector<double>> c;
  KernelState K;
};

DynamicSolution run_dynamic(const ChannelConfig& cfg,
                            const Grid1D& grid,
                            const KernelLibrary& kernels);

void write_dynamic_outputs(const ChannelConfig& cfg,
                           const DynamicSolution& sol,
                           const std::string& subdir = "");

} // namespace channel
