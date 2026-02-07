#pragma once

#include <channel/grid.hpp>
#include <channel/kernels.hpp>
#include <channel/types.hpp>

#include <optional>
#include <vector>

namespace channel {

struct ChargeSummary {
  double Q_vol = 0.0;      // ∫ ρ_free dz [C/m^2] in 1D per-area convention
  double Q_gate = 0.0;     // -ε0 ε(d) ψ'(d) [C/m^2]
  double Q_left = 0.0;     // +ε0 ε(0) ψ'(0) [C/m^2] (sign so Q_vol = Q_left + Q_gate)
  double maxwell_mismatch = 0.0; // |Q_vol - (Q_left + Q_gate)|
};

std::vector<double> compute_rho_free(const Grid1D& grid,
                                     const std::vector<Species>& species,
                                     const KernelState& K,
                                     const std::vector<std::vector<double>>& c,
                                     const std::vector<double>& alpha,
                                     const std::optional<RedoxParams>& redox);

ChargeSummary compute_charge(const Grid1D& grid,
                             const std::vector<Species>& species,
                             const KernelState& K,
                             const std::vector<std::vector<double>>& c,
                             const std::vector<double>& psi,
                             const std::vector<double>& alpha,
                             const std::optional<RedoxParams>& redox,
                             double psi0,
                             double psiD);

double alpha_bar_sites(const Grid1D& grid, const KernelState& K, const std::vector<double>& alpha);

} // namespace channel
