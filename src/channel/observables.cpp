#include <channel/observables.hpp>
#include <channel/constants.hpp>

#include <cmath>
#include <stdexcept>

namespace channel {

std::vector<double> compute_rho_free(const Grid1D& grid,
                                     const std::vector<Species>& species,
                                     const KernelState& K,
                                     const std::vector<std::vector<double>>& c,
                                     const std::vector<double>& alpha,
                                     const std::optional<RedoxParams>& redox) {
  const std::size_t nz = grid.n_cells();
  if (alpha.size() != nz) throw std::runtime_error("compute_rho_free: alpha length mismatch");
  if (c.size() != species.size()) throw std::runtime_error("compute_rho_free: c species mismatch");
  for (const auto& ci : c) if (ci.size() != nz) throw std::runtime_error("compute_rho_free: c[z] mismatch");

  std::vector<double> rho(nz, 0.0);
  for (std::size_t k = 0; k < nz; ++k) {
    double r = K.rho_base[k];
    for (std::size_t i = 0; i < species.size(); ++i) {
      r += species[i].q * c[i][k];
    }
    if (redox) {
      r += static_cast<double>(redox->sigmaP) * eCharge * K.ns[k] * alpha[k];
    }
    rho[k] = r;
  }
  return rho;
}

ChargeSummary compute_charge(const Grid1D& grid,
                             const std::vector<Species>& species,
                             const KernelState& K,
                             const std::vector<std::vector<double>>& c,
                             const std::vector<double>& psi,
                             const std::vector<double>& alpha,
                             const std::optional<RedoxParams>& redox,
                             double psi0,
                             double psiD) {
  const std::size_t nz = grid.n_cells();
  const double dz = grid.dz();

  if (psi.size() != nz) throw std::runtime_error("compute_charge: psi length mismatch");

  ChargeSummary out;

  auto rho = compute_rho_free(grid, species, K, c, alpha, redox);

  for (std::size_t k = 0; k < nz; ++k) out.Q_vol += rho[k] * dz;

  // Boundary derivatives at faces
  double dpsi0 = (psi[0] - psi0) / dz;
  double dpsiD = (psiD - psi[nz - 1]) / dz;

  out.Q_left = +epsilon0 * K.epsr[0] * dpsi0;
  out.Q_gate = -epsilon0 * K.epsr[nz - 1] * dpsiD;

  out.maxwell_mismatch = std::fabs(out.Q_vol - (out.Q_left + out.Q_gate));
  return out;
}

double alpha_bar_sites(const Grid1D& grid, const KernelState& K, const std::vector<double>& alpha) {
  const std::size_t nz = grid.n_cells();
  if (alpha.size() != nz) throw std::runtime_error("alpha_bar_sites: alpha length mismatch");
  const double dz = grid.dz();

  double num = 0.0;
  double den = 0.0;
  for (std::size_t k = 0; k < nz; ++k) {
    num += K.ns[k] * alpha[k] * dz;
    den += K.ns[k] * dz;
  }
  if (den <= 0.0) {
    // fallback: simple average
    double s = 0.0;
    for (double a : alpha) s += a;
    return s / static_cast<double>(alpha.size());
  }
  return num / den;
}

} // namespace channel
