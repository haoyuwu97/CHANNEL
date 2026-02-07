#include <channel/stationary.hpp>

#include <channel/constants.hpp>
#include <channel/functional.hpp>
#include <channel/io.hpp>
#include <channel/poisson.hpp>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <stdexcept>

namespace channel {

namespace fs = std::filesystem;

namespace {

double safe_exp(double x) {
  if (x > 700.0) x = 700.0;
  if (x < -700.0) x = -700.0;
  return std::exp(x);
}

double clamp01(double a) {
  if (a < 1e-15) return 1e-15;
  if (a > 1.0 - 1e-15) return 1.0 - 1e-15;
  return a;
}

std::size_t find_species(const std::vector<Species>& sp, const std::string& name) {
  for (std::size_t i = 0; i < sp.size(); ++i) if (sp[i].name == name) return i;
  throw std::runtime_error("Species not found: " + name);
}

// Compute c_i(z) from chemical potential equilibrium.
// For i != X: Boltzmann-like.
// For X with explicit coupling: solve Eq.(17) locally by Newton.
std::vector<std::vector<double>> update_concentrations(double T,
                                                       const Grid1D& grid,
                                                       const std::vector<Species>& species,
                                                       const KernelState& K,
                                                       const std::vector<double>& psi,
                                                       const std::vector<double>& alpha,
                                                       const std::optional<RedoxParams>& redox,
                                                       bool explicit_counterion_coupling,
                                                       const std::string& closure_mode) {
  const std::size_t nz = grid.n_cells();
  const double betaT = beta(T);
  const std::string mode = closure_mode;

  std::vector<std::vector<double>> c(species.size(), std::vector<double>(nz, 0.0));

  int idxX = -1;
  if (redox) idxX = static_cast<int>(find_species(species, redox->counterion));

  for (std::size_t i = 0; i < species.size(); ++i) {
    for (std::size_t k = 0; k < nz; ++k) {
      double ci = 0.0;

      if (species[i].c_res > 0.0) {
        if (mode == "c") {
          // Mode C: c = c_res exp[-(phi_ex + β q ψ)]
          double logc = std::log(species[i].c_res)
                      - (K.sp[i].phi_ex[k] + betaT * species[i].q * psi[k]);
          if (logc > 700.0) logc = 700.0;
          if (logc < -700.0) logc = -700.0;
          ci = std::exp(logc);
        } else {
          // Ω-based (Mode A/B): c = h c_res exp[-β(Uex + q ψ)]
          double h = K.sp[i].hi[k];
          if (h <= 0.0) throw std::runtime_error("update_concentrations: h_i <= 0");
          double expo = -betaT * (K.sp[i].Uex[k] + species[i].q * psi[k]);
          double logc = std::log(h) + std::log(species[i].c_res) + expo;
          if (logc > 700.0) logc = 700.0;
          if (logc < -700.0) logc = -700.0;
          ci = std::exp(logc);
        }
      }

      if (ci < 0.0) ci = 0.0;
      c[i][k] = ci;
    }
  }

  // Counterion coupling term from ω_redox if requested.
  if (redox && explicit_counterion_coupling && idxX >= 0) {
    for (std::size_t k = 0; k < nz; ++k) {
      double ns = K.ns[k];
      double a = alpha[k];
      double h = K.sp[idxX].hi[k];
      double cres = species[idxX].c_res;
      double U = K.sp[idxX].Uex[k] + species[idxX].q * psi[k];

      if (ns <= 0.0 || a <= 0.0) continue;

      // Solve g(c)= ln(c/(h cres)) + β U - ns a / c = 0.
      double c0 = c[idxX][k];
      if (c0 < 1e-30) c0 = h * cres; // fallback
      double ccur = c0;

      for (int it = 0; it < 50; ++it) {
        if (ccur < 1e-30) ccur = 1e-30;
        double g = std::log(ccur / (h * cres)) + betaT * U - ns * a / ccur;
        double gp = (1.0 / ccur) + (ns * a) / (ccur * ccur);

        double step = g / gp;
        // Damped Newton to keep positivity
        double cnew = ccur - step;
        if (cnew <= 0.0) cnew = 0.5 * ccur;

        if (std::fabs(g) < 1e-10) {
          ccur = cnew;
          break;
        }
        ccur = cnew;
      }

      c[idxX][k] = ccur;
    }
  }

  return c;
}

std::vector<double> update_alpha_no_feedback(double T,
                                             const Grid1D& grid,
                                             const KernelState& K,
                                             const std::vector<Species>& species,
                                             const std::vector<double>& psi,
                                             const std::vector<std::vector<double>>& c,
                                             const RedoxParams& redox) {
  const std::size_t nz = grid.n_cells();
  const double betaT = beta(T);

  std::size_t idxX = find_species(species, redox.counterion);

  std::vector<double> a(nz, 0.0);
  for (std::size_t k = 0; k < nz; ++k) {
    double cX = c[idxX][k];
    if (cX < 1e-300) cX = 1e-300;

    double log_ratio = std::log(redox.KX) + std::log(cX)
                     - betaT * (redox.deltaG0 + static_cast<double>(redox.sigmaP) * eCharge * psi[k]);
    double ak;
    if (log_ratio > 50.0) {
      ak = 1.0;
    } else if (log_ratio < -50.0) {
      ak = 0.0;
    } else {
      ak = 1.0 / (1.0 + std::exp(-log_ratio));
    }
    a[k] = clamp01(ak);
  }
  return a;
}

std::vector<double> update_alpha_with_feedback(double T,
                                               const Grid1D& grid,
                                               const KernelState& K,
                                               const std::vector<Species>& species,
                                               const std::vector<double>& psi,
                                               const std::vector<std::vector<double>>& c,
                                               const std::vector<double>& alpha_old,
                                               const RedoxParams& redox) {
  const std::size_t nz = grid.n_cells();
  const double dz = grid.dz();
  const double betaT = beta(T);
  const double kBT = kBoltzmann * T;

  std::size_t idxX = find_species(species, redox.counterion);

  // Precompute |ψ'|^2 at cell centers
  std::vector<double> dpsi2(nz, 0.0);
  for (std::size_t k = 0; k < nz; ++k) {
    double dpsi = 0.0;
    if (k == 0) dpsi = (psi[1] - psi[0]) / dz;
    else if (k + 1 == nz) dpsi = (psi[nz - 1] - psi[nz - 2]) / dz;
    else dpsi = (psi[k + 1] - psi[k - 1]) / (2.0 * dz);
    dpsi2[k] = dpsi * dpsi;
  }

  std::vector<double> a = alpha_old;

  for (std::size_t k = 0; k < nz; ++k) {
    double ns = K.ns[k];
    if (ns <= 0.0) {
      a[k] = 0.0;
      continue;
    }

    double cX = c[idxX][k];
    if (cX < 1e-300) cX = 1e-300;

    // Feedback scalar F = (1/ns) Σ_i [c_i dUex_i/dα - kBT c_i dlnh_i/dα] - (ε0/(2ns)) depsr/dα |ψ'|^2
    double F = 0.0;
    for (std::size_t i = 0; i < species.size(); ++i) {
      F += c[i][k] * K.sp[i].dUex_dalpha[k] - kBT * c[i][k] * K.sp[i].dlnhi_dalpha[k];
    }
    F /= ns;
    F += -(epsilon0 / (2.0 * ns)) * K.depsr_dalpha[k] * dpsi2[k];
    // Mode B: include ω_extra contribution through its local derivative.
    if (!K.domega_extra_dalpha.empty()) {
      F += K.domega_extra_dalpha[k] / ns;
    }

    // Solve f(α) = kBT ln(α/(1-α)) + ΔG0 - kBT ln(KX cX) + σP e ψ + F = 0
    // Newton in α with bounds.
    double alpha_k = clamp01(a[k]);

    const double base = redox.deltaG0 - kBT * std::log(redox.KX * cX)
                        + static_cast<double>(redox.sigmaP) * eCharge * psi[k]
                        + F;

    for (int it = 0; it < 50; ++it) {
      alpha_k = clamp01(alpha_k);

      double f = kBT * std::log(alpha_k / (1.0 - alpha_k)) + base;
      double fp = kBT * (1.0 / alpha_k + 1.0 / (1.0 - alpha_k));

      double step = f / fp;
      double anew = alpha_k - step;

      // clamp and damp
      anew = clamp01(anew);
      if (std::fabs(f) < 1e-10) {
        alpha_k = anew;
        break;
      }
      alpha_k = anew;
    }
    a[k] = alpha_k;
  }

  return a;
}

double max_abs_diff(const std::vector<double>& a, const std::vector<double>& b) {
  if (a.size() != b.size()) throw std::runtime_error("max_abs_diff: size mismatch");
  double m = 0.0;
  for (std::size_t i = 0; i < a.size(); ++i) m = std::max(m, std::fabs(a[i] - b[i]));
  return m;
}

} // namespace

StationarySolution solve_stationary(const ChannelConfig& cfg,
                                    const Grid1D& grid,
                                    const KernelLibrary& kernels,
                                    double VG_override,
                                    bool compute_capacitance) {
  StationarySolution sol;
  sol.VG = VG_override;

  const std::size_t nz = grid.n_cells();

  // Initial guess
  sol.psi.assign(nz, 0.0);
  for (std::size_t k = 0; k < nz; ++k) {
    double z = grid.zc()[k];
    sol.psi[k] = (z / grid.length()) * sol.VG;
  }

  sol.alpha.assign(nz, 0.0);
  if (cfg.redox) {
    for (std::size_t k = 0; k < nz; ++k) sol.alpha[k] = clamp01(cfg.alpha_init);
  }

  // Picard iterations
  std::vector<double> psi_old = sol.psi;
  std::vector<double> alpha_old = sol.alpha;

  for (int it = 0; it < cfg.max_iter; ++it) {
    sol.iterations = it + 1;

    // Evaluate kernels for current α(z)
    sol.K = evaluate_kernels(kernels, grid, cfg.species, sol.alpha);

    // Update ion concentrations
    sol.c = update_concentrations(cfg.T, grid, cfg.species, sol.K, sol.psi, sol.alpha, cfg.redox,
                                  cfg.explicit_counterion_coupling, cfg.closure_mode);

    // Update ψ from Poisson
    auto rho = compute_rho_free(grid, cfg.species, sol.K, sol.c, sol.alpha, cfg.redox);
    auto psi_new = solve_poisson_dirichlet(grid, sol.K.epsr, rho, 0.0, sol.VG);

    // Mix ψ
    for (std::size_t k = 0; k < nz; ++k) {
      sol.psi[k] = (1.0 - cfg.damping) * sol.psi[k] + cfg.damping * psi_new[k];
    }

    // Update α if redox
    if (cfg.redox) {
      std::vector<double> alpha_new;
      if (cfg.enable_feedback) {
        alpha_new = update_alpha_with_feedback(cfg.T, grid, sol.K, cfg.species, sol.psi, sol.c, sol.alpha, *cfg.redox);
      } else {
        alpha_new = update_alpha_no_feedback(cfg.T, grid, sol.K, cfg.species, sol.psi, sol.c, *cfg.redox);
      }
      for (std::size_t k = 0; k < nz; ++k) {
        sol.alpha[k] = (1.0 - cfg.damping) * sol.alpha[k] + cfg.damping * alpha_new[k];
        sol.alpha[k] = clamp01(sol.alpha[k]);
      }
    }

    // Convergence
    double dpsi = max_abs_diff(sol.psi, psi_old);
    double dalpha = max_abs_diff(sol.alpha, alpha_old);
    psi_old = sol.psi;
    alpha_old = sol.alpha;

    if (dpsi < cfg.tol && dalpha < cfg.tol) {
      sol.converged = true;
      break;
    }
  }

  // Final recompute for consistency (kernels, c)
  sol.K = evaluate_kernels(kernels, grid, cfg.species, sol.alpha);
  sol.c = update_concentrations(cfg.T, grid, cfg.species, sol.K, sol.psi, sol.alpha, cfg.redox,
                                cfg.explicit_counterion_coupling, cfg.closure_mode);
  sol.charge = compute_charge(grid, cfg.species, sol.K, sol.c, sol.psi, sol.alpha, cfg.redox, 0.0, sol.VG);
  sol.Omega = compute_Omega(cfg.T, grid, cfg.species, sol.K, sol.c, sol.psi, sol.alpha, cfg.redox,
                            cfg.closure_mode).Omega;

  // Capacitance estimate via finite difference
  if (compute_capacitance && cfg.compute_capacitance) {
    double dV = cfg.cap_deltaV;
    StationarySolution plus = solve_stationary(cfg, grid, kernels, sol.VG + dV, false);
    StationarySolution minus = solve_stationary(cfg, grid, kernels, sol.VG - dV, false);
    sol.has_capacitance = true;
    sol.C_eq = (plus.charge.Q_gate - minus.charge.Q_gate) / (2.0 * dV);
  }

  return sol;
}

void write_stationary_outputs(const ChannelConfig& cfg,
                              const Grid1D& grid,
                              const std::vector<Species>& species,
                              const StationarySolution& sol,
                              const std::string& subdir) {
  fs::path outdir(cfg.output_dir);
  if (!subdir.empty()) outdir /= subdir;
  ensure_dir(outdir.string());

  const std::vector<double>& z = grid.zc();
  const std::size_t nz = grid.n_cells();

  // ψ(z)
  {
    std::vector<std::vector<double>> cols = {z, sol.psi};
    write_table((outdir / "psi.dat").string(), {"z", "psi"}, cols, "z [m], psi [V]");
  }

  // α(z) if present
  if (cfg.redox) {
    std::vector<std::vector<double>> cols = {z, sol.alpha};
    write_table((outdir / "alpha.dat").string(), {"z", "alpha"}, cols, "z [m], alpha [-]");
  }

  // c_i(z)
  for (std::size_t i = 0; i < species.size(); ++i) {
    std::vector<std::vector<double>> cols = {z, sol.c[i]};
    write_table((outdir / ("c_" + species[i].name + ".dat")).string(),
                {"z", "c"}, cols, "z [m], c [1/m^3]");
  }

  // Kernel snippets
  {
    std::vector<std::vector<double>> cols = {z, sol.K.epsr, sol.K.ns, sol.K.rho_base};
    write_table((outdir / "kernels_core.dat").string(),
                {"z", "epsr", "ns", "rho_base"}, cols,
                "epsr [-], ns [1/m^3], rho_base [C/m^3]");
  }

  // summary.json
  ResultsIndex idx;
  idx.mode = "stationary";
  idx.config_used = "config_used.ini";
  idx.summary["VG"] = std::to_string(sol.VG);
  idx.summary["converged"] = sol.converged ? "true" : "false";
  idx.summary["iterations"] = std::to_string(sol.iterations);
  idx.summary["Q_gate_C_per_m2"] = std::to_string(sol.charge.Q_gate);
  idx.summary["Q_vol_C_per_m2"] = std::to_string(sol.charge.Q_vol);
  idx.summary["maxwell_mismatch_C_per_m2"] = std::to_string(sol.charge.maxwell_mismatch);
  idx.summary["Omega_J_per_m2"] = std::to_string(sol.Omega);
  if (sol.has_capacitance) idx.summary["C_eq_F_per_m2"] = std::to_string(sol.C_eq);

  idx.datasets["psi"] = DatasetMeta{"psi.dat", {"z", "psi"}, "Electrostatic potential"};
  if (cfg.redox) idx.datasets["alpha"] = DatasetMeta{"alpha.dat", {"z", "alpha"}, "Oxidation fraction"};
  for (const auto& sp : species) {
    idx.datasets["c_" + sp.name] = DatasetMeta{"c_" + sp.name + ".dat", {"z", "c"}, "Concentration of " + sp.name};
  }
  idx.datasets["kernels_core"] = DatasetMeta{"kernels_core.dat", {"z", "epsr", "ns", "rho_base"}, "Core kernel profiles"};

  write_results_json(outdir.string(), idx);
}

} // namespace channel
