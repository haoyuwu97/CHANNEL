#include <channel/functional.hpp>
#include <channel/constants.hpp>

#include <cmath>
#include <limits>
#include <stdexcept>

namespace channel {

FunctionalContrib compute_Omega(double T,
                                const Grid1D& grid,
                                const std::vector<Species>& species,
                                const KernelState& K,
                                const std::vector<std::vector<double>>& c,
                                const std::vector<double>& psi,
                                const std::vector<double>& alpha,
                                const std::optional<RedoxParams>& redox,
                                const std::string& closure_mode) {
  const std::size_t nz = grid.n_cells();
  const double dz = grid.dz();

  if (psi.size() != nz || alpha.size() != nz) {
    throw std::runtime_error("compute_Omega: psi/alpha length mismatch with grid");
  }
  if (c.size() != species.size()) {
    throw std::runtime_error("compute_Omega: c species dimension mismatch");
  }
  for (const auto& ci : c) {
    if (ci.size() != nz) throw std::runtime_error("compute_Omega: c[z] length mismatch");
  }

  const std::string mode = closure_mode;
  if (mode == "c") {
    // μ-closure: Ω is not defined by default (audits that require Ω should be skipped).
    const double nan = std::numeric_limits<double>::quiet_NaN();
    FunctionalContrib out;
    out.Omega = nan;
    out.Omega_field = nan;
    out.Omega_ion = nan;
    out.Omega_redox = nan;
    out.Omega_extra = nan;
    return out;
  }

  const double kBT = kBoltzmann * T;

  FunctionalContrib out;

  // Identify counter-ion index for ω_redox
  int idxX = -1;
  if (redox) {
    for (std::size_t i = 0; i < species.size(); ++i) {
      if (species[i].name == redox->counterion) {
        idxX = static_cast<int>(i);
        break;
      }
    }
    if (idxX < 0) {
      throw std::runtime_error("compute_Omega: redox enabled but counterion not found in species list");
    }
  }

  for (std::size_t k = 0; k < nz; ++k) {
    // ψ' via finite difference (approx)
    double dpsi = 0.0;
    if (k == 0) {
      dpsi = (psi[1] - psi[0]) / dz;
    } else if (k + 1 == nz) {
      dpsi = (psi[nz - 1] - psi[nz - 2]) / dz;
    } else {
      dpsi = (psi[k + 1] - psi[k - 1]) / (2.0 * dz);
    }

    double rho = K.rho_base[k];
    for (std::size_t i = 0; i < species.size(); ++i) {
      rho += species[i].q * c[i][k];
    }
    if (redox) {
      rho += static_cast<double>(redox->sigmaP) * eCharge * K.ns[k] * alpha[k];
    }

    double omega_field = -0.5 * epsilon0 * K.epsr[k] * dpsi * dpsi + psi[k] * rho;

    double omega_ion = 0.0;
    for (std::size_t i = 0; i < species.size(); ++i) {
      double ci = c[i][k];
      if (ci < 1e-300) ci = 1e-300;
      double hi = K.sp[i].hi[k];
      if (hi < 1e-300) hi = 1e-300;
      double cres = species[i].c_res;
      if (cres < 1e-300) cres = 1e-300;

      omega_ion += kBT * ci * (std::log(ci / (hi * cres)) - 1.0) + ci * K.sp[i].Uex[k];
    }

    double omega_redox = 0.0;
    if (redox) {
      double a = alpha[k];
      if (a < 1e-15) a = 1e-15;
      if (a > 1.0 - 1e-15) a = 1.0 - 1e-15;
      double cX = c[idxX][k];
      if (cX < 1e-300) cX = 1e-300;

      omega_redox = K.ns[k] * (kBT * (a * std::log(a) + (1.0 - a) * std::log(1.0 - a))
                              + a * redox->deltaG0
                              - kBT * a * std::log(redox->KX * cX));
    }

    out.Omega_field += omega_field * dz;
    out.Omega_ion += omega_ion * dz;
    out.Omega_redox += omega_redox * dz;

    // Mode B: ω_extra(z;α) is included as an additive energy density.
    // For mode A, omega_extra is typically identically zero.
    if (!K.omega_extra.empty()) {
      out.Omega_extra += K.omega_extra[k] * dz;
    }
  }

  out.Omega = out.Omega_field + out.Omega_ion + out.Omega_redox + out.Omega_extra;
  return out;
}

} // namespace channel
