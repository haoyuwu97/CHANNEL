#include <channel/verify.hpp>

#include <cmath>
#include <limits>
#include <stdexcept>

namespace channel {

VerificationReport verify_stationary(const ChannelConfig& cfg,
                                     const Grid1D& grid,
                                     const KernelLibrary& kernels,
                                     const StationarySolution& sol) {
  VerificationReport r;
  r.maxwell_mismatch = sol.charge.maxwell_mismatch;
  r.Q_gate = sol.charge.Q_gate;

  // Maxwell (Gauss-law) consistency
  const double tol_Q = 1e-6 * std::max(1.0, std::fabs(sol.charge.Q_gate));
  r.maxwell_ok = (r.maxwell_mismatch <= tol_Q);

  // μ-closure (Mode C): Ω is not defined by default, so Maxwell/energy audits are not applicable.
  if (cfg.closure_mode == "c") {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    r.dOmega_dVG = nan;
    r.energy_rel_error = nan;
    r.energy_ok = true; // treat as N/A
    return r;
  }

  // Energy Maxwell relation: dΩ/dVG = Q (per-area convention).
  const double dV = (cfg.cap_deltaV > 0.0) ? cfg.cap_deltaV : 1e-3;

  StationarySolution plus = solve_stationary(cfg, grid, kernels, sol.VG + dV, false);
  StationarySolution minus = solve_stationary(cfg, grid, kernels, sol.VG - dV, false);

  r.dOmega_dVG = (plus.Omega - minus.Omega) / (2.0 * dV);

  double denom = std::max(1.0, std::fabs(r.Q_gate));
  r.energy_rel_error = std::fabs(r.dOmega_dVG - r.Q_gate) / denom;

  r.energy_ok = (r.energy_rel_error < 1e-2); // loose by default (discretization effects)

  return r;
}

} // namespace channel
