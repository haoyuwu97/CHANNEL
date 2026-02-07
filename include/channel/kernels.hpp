#pragma once

#include <channel/config.hpp>
#include <channel/grid.hpp>
#include <channel/profile.hpp>
#include <channel/types.hpp>

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace channel {

struct KernelField {
  std::optional<Profile1D> profile;
  std::optional<FieldZAlpha> field;

  static KernelField constant(double value, const Grid1D& grid);
  static KernelField from_profile(Profile1D p);
  static KernelField from_field(FieldZAlpha f);

  bool is_alpha_dependent() const { return field.has_value(); }

  // Evaluate on z_query using the local alpha_field(z_query).
  // If dout_dalpha is not null, returns ∂f/∂α (local derivative) on z_query.
  void eval(const std::vector<double>& z_query,
            const std::vector<double>& alpha_field,
            std::vector<double>& out,
            std::vector<double>* dout_dalpha = nullptr) const;
};

struct SpeciesKernels {
  std::string name;
  KernelField hi;          // dimensionless
  KernelField delta_mu0;   // [J]

  // Mode C: dimensionless excess potential ϕ_ex(z;α) used to close μ without an explicit Ω.
  // Convention: μ_i/(kBT) = ln(c_i/c_i^{res}) + ϕ_ex,i(z;α) + β q_i ψ.
  // Therefore, equilibrium with a reservoir reference implies
  //   c_i(z) = c_i^{res} exp[-(ϕ_ex,i + β q_i ψ)].
  KernelField phi_ex;
};

struct KernelLibrary {
  double T = 298.15;
  double epsr_res = 78.0;
  std::string parameterization = "a"; // "a" or "b" (lowercase)

  // Closure mode: a/b => Ω-based, c => μ-based.
  std::string closure_mode = "a";

  KernelField epsr;
  KernelField ns;
  KernelField rho_base;

  // Mode B: additional (unknown / surrogate) energy density term ω_extra(z;α) [J/m^3].
  KernelField omega_extra;

  std::vector<SpeciesKernels> sp;

  std::size_t species_index(const std::string& name) const;
};

struct KernelState {
  std::vector<double> epsr;
  std::vector<double> depsr_dalpha;

  std::vector<double> ns;
  std::vector<double> rho_base;

  // Mode B energy term ω_extra and derivative.
  std::vector<double> omega_extra;
  std::vector<double> domega_extra_dalpha;

  struct SpState {
    std::vector<double> hi;
    std::vector<double> dlnhi_dalpha;

    std::vector<double> delta_mu0;
    std::vector<double> ddelta_mu0_dalpha;

    std::vector<double> Uex;
    std::vector<double> dUex_dalpha;

    // Mode C potentials (dimensionless) and derivatives.
    std::vector<double> phi_ex;
    std::vector<double> dphi_ex_dalpha;
  };
  std::vector<SpState> sp;
};

KernelLibrary build_kernel_library(const ChannelConfig& cfg, const Grid1D& grid);

KernelState evaluate_kernels(const KernelLibrary& lib,
                             const Grid1D& grid,
                             const std::vector<Species>& species,
                             const std::vector<double>& alpha_field);

} // namespace channel
