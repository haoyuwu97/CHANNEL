#pragma once

#include <channel/grid.hpp>
#include <channel/kernels.hpp>
#include <channel/types.hpp>

#include <optional>
#include <vector>

namespace channel {

struct FunctionalContrib {
  double Omega = 0.0;
  double Omega_field = 0.0;
  double Omega_ion = 0.0;
  double Omega_redox = 0.0;
  double Omega_extra = 0.0;
};

// Compute the grand potential Ω for the current fields, in the discretized 1D geometry.
// This is mainly used for verification (Maxwell relation, energy decay).
FunctionalContrib compute_Omega(double T,
                                const Grid1D& grid,
                                const std::vector<Species>& species,
                                const KernelState& K,
                                const std::vector<std::vector<double>>& c, // [nspecies][nz]
                                const std::vector<double>& psi,            // [nz]
                                const std::vector<double>& alpha,          // [nz]
                                const std::optional<RedoxParams>& redox,
                                const std::string& closure_mode);

} // namespace channel
