#include <channel/dynamic.hpp>

#include <channel/device.hpp>

#include <channel/constants.hpp>
#include <channel/functional.hpp>
#include <channel/io.hpp>
#include <channel/observables.hpp>
#include <channel/poisson.hpp>
#include <channel/stationary.hpp>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>
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

double bernoulli(double x) {
  // B(x) = x/(exp(x)-1) with overflow-safe asymptotics.
  const double ax = std::fabs(x);
  if (ax < 1e-6) {
    double x2 = x * x;
    return 1.0 - 0.5 * x + x2 / 12.0 - x2 * x2 / 720.0;
  }
  if (x > 50.0) {
    // exp(x) dominates
    return x * std::exp(-x);
  }
  if (x < -50.0) {
    // exp(x) ~ 0, so exp(x)-1 ~ -1
    return -x;
  }
  return x / std::expm1(x);
}

static void thomas(std::vector<double>& a, std::vector<double>& b, std::vector<double>& c, std::vector<double>& d) {
  const std::size_t n = b.size();
  for (std::size_t i = 1; i < n; ++i) {
    double m = a[i] / b[i - 1];
    b[i] -= m * c[i - 1];
    d[i] -= m * d[i - 1];
  }
  d[n - 1] /= b[n - 1];
  for (std::size_t i = n - 1; i-- > 0;) {
    d[i] = (d[i] - c[i] * d[i + 1]) / b[i];
  }
}

struct Waveform {
  std::string kind;
  double VG0 = 0.0;
  double VG1 = 0.0;
  double freq = 0.0;
  std::vector<double> t;
  std::vector<double> VG;

  double eval(double time) const {
    if (kind == "step") {
      return (time <= 0.0) ? VG0 : VG1;
    }
    if (kind == "sine") {
      return VG0 + VG1 * std::sin(2.0 * pi * freq * time);
    }
    if (kind == "triangle") {
      double phase = std::fmod(freq * time, 1.0);
      if (phase < 0.0) phase += 1.0;
      // triangle in [-1,1]
      double tri = (phase < 0.5) ? (4.0 * phase - 1.0) : (-4.0 * phase + 3.0);
      return VG0 + VG1 * tri;
    }
    if (kind == "file") {
      if (t.empty()) return VG0;
      if (time <= t.front()) return VG.front();
      if (time >= t.back()) return VG.back();
      auto it = std::upper_bound(t.begin(), t.end(), time);
      std::size_t i = static_cast<std::size_t>(it - t.begin());
      std::size_t i0 = i - 1, i1 = i;
      double t0 = t[i0], t1 = t[i1];
      double s = (time - t0) / (t1 - t0);
      return (1.0 - s) * VG[i0] + s * VG[i1];
    }
    return VG0;
  }
};

Waveform build_waveform(const ChannelConfig& cfg) {
  Waveform w;
  w.kind = cfg.waveform;
  w.VG0 = cfg.VG0;
  w.VG1 = cfg.VG1;
  w.freq = cfg.freq;
  if (w.kind == "file") {
    std::ifstream f(cfg.waveform_file);
    if (!f) throw std::runtime_error("Cannot open waveform_file: " + cfg.waveform_file);
    std::string line;
    while (std::getline(f, line)) {
      if (line.empty()) continue;
      if (line[0] == '#' || line[0] == ';') continue;
      std::stringstream ss(line);
      double t = 0.0, v = 0.0;
      if (!(ss >> t >> v)) continue;
      w.t.push_back(t);
      w.VG.push_back(v);
    }
    if (w.t.size() < 2) throw std::runtime_error("waveform_file needs at least 2 points");
  }
  return w;
}

std::size_t find_species(const std::vector<Species>& sp, const std::string& name) {
  for (std::size_t i = 0; i < sp.size(); ++i) if (sp[i].name == name) return i;
  throw std::runtime_error("Species not found: " + name);
}

} // namespace

DynamicSolution run_dynamic(const ChannelConfig& cfg,
                            const Grid1D& grid,
                            const KernelLibrary& kernels) {
  if (!cfg.dynamic_enabled) {
    throw std::runtime_error("run_dynamic called but dynamic.enabled is false");
  }

  const std::size_t nz = grid.n_cells();
  const double dt = cfg.dt;
  const int n_steps = static_cast<int>(std::ceil(cfg.t_end / dt));
  Waveform wf = build_waveform(cfg);

  // Initial condition: equilibrium at VG0 (stationary solve)
  StationarySolution init = solve_stationary(cfg, grid, kernels, cfg.VG0, false);

  DynamicSolution sol;
  sol.psi = init.psi;
  sol.alpha = init.alpha;
  sol.c = init.c;
  sol.K = init.K;

  sol.t.reserve(n_steps + 1);
  sol.VG.reserve(n_steps + 1);
  sol.Q_gate.reserve(n_steps + 1);
  sol.Q_vol.reserve(n_steps + 1);
  sol.alpha_bar.reserve(n_steps + 1);
  sol.Omega.reserve(n_steps + 1);

  // Record at t=0
  double t = 0.0;
  double VG = wf.eval(t);
  // Ensure ψ consistent with VG at t=0 (init uses VG0, but waveform might define differently)
  {
    sol.K = evaluate_kernels(kernels, grid, cfg.species, sol.alpha);
    auto rho = compute_rho_free(grid, cfg.species, sol.K, sol.c, sol.alpha, cfg.redox);
    sol.psi = solve_poisson_dirichlet(grid, sol.K.epsr, rho, 0.0, VG);
  }
  auto charge0 = compute_charge(grid, cfg.species, sol.K, sol.c, sol.psi, sol.alpha, cfg.redox, 0.0, VG);
  double Omega0 = compute_Omega(cfg.T, grid, cfg.species, sol.K, sol.c, sol.psi, sol.alpha, cfg.redox,
                                cfg.closure_mode).Omega;

  sol.t.push_back(t);
  sol.VG.push_back(VG);
  sol.Q_gate.push_back(charge0.Q_gate);
  sol.Q_vol.push_back(charge0.Q_vol);
  sol.alpha_bar.push_back(alpha_bar_sites(grid, sol.K, sol.alpha));
  sol.Omega.push_back(Omega0);

  // Helper: reaction update
  auto do_reaction = [&](const std::vector<double>& alpha_in,
                         const std::vector<std::vector<double>>& c_in,
                         const std::vector<double>& psi_in) -> std::vector<double> {
    if (!cfg.redox) return std::vector<double>(nz, 0.0);
    const RedoxParams& r = *cfg.redox;
    std::size_t idxX = find_species(cfg.species, r.counterion);
    const double betaT = beta(cfg.T);
    const double kBT = kBoltzmann * cfg.T;

    // Kernels at the current α-field (needed for feedback terms).
    KernelState Kreac = evaluate_kernels(kernels, grid, cfg.species, alpha_in);

    // Precompute |ψ'|^2 for feedback if needed.
    std::vector<double> dpsi2(nz, 0.0);
    if (cfg.enable_feedback) {
      const double dz = grid.dz();
      for (std::size_t k = 0; k < nz; ++k) {
        double dpsi = 0.0;
        if (k == 0) dpsi = (psi_in[1] - psi_in[0]) / dz;
        else if (k + 1 == nz) dpsi = (psi_in[nz - 1] - psi_in[nz - 2]) / dz;
        else dpsi = (psi_in[k + 1] - psi_in[k - 1]) / (2.0 * dz);
        dpsi2[k] = dpsi * dpsi;
      }
    }

    std::vector<double> out(alpha_in);
    for (std::size_t k = 0; k < nz; ++k) {
      double cX = c_in[idxX][k];
      if (cX < 1e-300) cX = 1e-300;

      // Optional α-feedback term (Ω-based closure only).
      double F = 0.0; // [J] per site
      if (cfg.enable_feedback) {
        double ns = Kreac.ns[k];
        if (ns > 0.0) {
          for (std::size_t i = 0; i < cfg.species.size(); ++i) {
            F += c_in[i][k] * Kreac.sp[i].dUex_dalpha[k] - kBT * c_in[i][k] * Kreac.sp[i].dlnhi_dalpha[k];
          }
          F /= ns;
          F += -(epsilon0 / (2.0 * ns)) * Kreac.depsr_dalpha[k] * dpsi2[k];
          // Mode B: include ω_extra contribution through its local derivative.
          if (!Kreac.domega_extra_dalpha.empty()) {
            F += Kreac.domega_extra_dalpha[k] / ns;
          }
        }
      }

      double log_ratio = std::log(r.KX) + std::log(cX)
                       - betaT * (r.deltaG0 + static_cast<double>(r.sigmaP) * eCharge * psi_in[k] + F);
      double alpha_eq;
      if (log_ratio > 50.0) alpha_eq = 1.0;
      else if (log_ratio < -50.0) alpha_eq = 0.0;
      else alpha_eq = 1.0 / (1.0 + std::exp(-log_ratio));

      if (r.koff <= 0.0) {
        out[k] = clamp01(alpha_eq);
      } else {
        double kon;
        {
          // detailed-balance: kon/koff = KX exp[-β(ΔG0 + σP e ψ + F)]
          double log_kon = std::log(r.koff) + std::log(r.KX)
                           - betaT * (r.deltaG0 + static_cast<double>(r.sigmaP) * eCharge * psi_in[k] + F);
          if (log_kon > 700.0) log_kon = 700.0;
          if (log_kon < -700.0) log_kon = -700.0;
          kon = std::exp(log_kon);
        }
        double lambda = kon * cX + r.koff;
        if (lambda <= 0.0) {
          out[k] = clamp01(alpha_eq);
        } else {
          double alpha_inf = (kon * cX) / lambda;
          out[k] = clamp01(alpha_inf + (alpha_in[k] - alpha_inf) * std::exp(-lambda * dt));
        }
      }
    }
    return out;
  };

  // NP transport update per species (implicit SG)
  auto transport_step = [&](const KernelState& Kstep,
                            const std::vector<double>& psi_step,
                            const std::vector<double>& alpha_old,
                            const std::vector<double>& alpha_new,
                            const std::vector<std::vector<double>>& c_old_in) -> std::vector<std::vector<double>> {
    std::vector<std::vector<double>> c_new = c_old_in;

    const double betaT = beta(cfg.T);
    const std::string mode = cfg.closure_mode;

    int idxX = -1;
    if (cfg.redox) idxX = static_cast<int>(find_species(cfg.species, cfg.redox->counterion));

    for (std::size_t i = 0; i < cfg.species.size(); ++i) {
      const Species& sp = cfg.species[i];
      if (!sp.mobile || sp.D <= 0.0) {
        continue;
      }

      // Source term S_i
      std::vector<double> S(nz, 0.0);
      if (cfg.redox && (static_cast<int>(i) == idxX) && cfg.reaction_first) {
        for (std::size_t k = 0; k < nz; ++k) {
          S[k] = -Kstep.ns[k] * (alpha_new[k] - alpha_old[k]) / dt;
        }
      }

      // Drift potential Φ (dimensionless):
      //  - Ω-based (A/B): Φ = β(Uex + qψ) - ln h
      //  - μ-closure (C): Φ = ϕ_ex + β qψ
      std::vector<double> Phi(nz, 0.0);
      for (std::size_t k = 0; k < nz; ++k) {
        if (mode == "c") {
          Phi[k] = Kstep.sp[i].phi_ex[k] + betaT * sp.q * psi_step[k];
        } else {
          double h = Kstep.sp[i].hi[k];
          if (h <= 0.0) h = 1e-300;
          Phi[k] = betaT * (Kstep.sp[i].Uex[k] + sp.q * psi_step[k]) - std::log(h);
        }
      }

      // Face coefficients between cells j and j+1
      std::vector<double> A(nz - 1, 0.0);
      std::vector<double> B(nz - 1, 0.0);
      for (std::size_t j = 0; j + 1 < nz; ++j) {
        double dPhi = Phi[j + 1] - Phi[j];
        A[j] = sp.D / grid.dz() * bernoulli(dPhi);
        B[j] = sp.D / grid.dz() * bernoulli(-dPhi);
      }

      // Left boundary Dirichlet c(0)=c_res, with distance dz/2.
      // Reservoir reference: Φ_L = 0 by convention (Uex_res=0, h_res=1, ψ(0)=0).
      double Phi_L = 0.0;
      double dPhiL = Phi[0] - Phi_L;
      double A_L = sp.D / (0.5 * grid.dz()) * bernoulli(dPhiL);
      double B_L = sp.D / (0.5 * grid.dz()) * bernoulli(-dPhiL);
      double c_L = sp.c_res;

      std::vector<double> a(nz, 0.0);
      std::vector<double> b(nz, 0.0);
      std::vector<double> ccoef(nz, 0.0);
      std::vector<double> rhs(nz, 0.0);

      for (std::size_t k = 0; k < nz; ++k) {
        rhs[k] = c_old_in[i][k] + dt * S[k];
        b[k] = 1.0;
      }

      // k=0
      {
        double B_right = B[0];
        double A_right = A[0];
        b[0] = 1.0 + (dt / grid.dz()) * (B_right + A_L);
        ccoef[0] = -(dt / grid.dz()) * A_right;
        rhs[0] += (dt / grid.dz()) * B_L * c_L;
      }

      // interior
      for (std::size_t k = 1; k + 1 < nz; ++k) {
        double B_right = B[k];
        double A_left  = A[k - 1];
        b[k] = 1.0 + (dt / grid.dz()) * (B_right + A_left);
        a[k] = -(dt / grid.dz()) * B[k - 1];
        ccoef[k] = -(dt / grid.dz()) * A[k];
      }

      // last
      {
        std::size_t k = nz - 1;
        double A_left = A[k - 1];
        b[k] = 1.0 + (dt / grid.dz()) * (A_left);
        a[k] = -(dt / grid.dz()) * B[k - 1];
        ccoef[k] = 0.0; // no flux right
      }

      thomas(a, b, ccoef, rhs);
      for (std::size_t k = 0; k < nz; ++k) {
        if (rhs[k] < 0.0) rhs[k] = 0.0;
        c_new[i][k] = rhs[k];
      }
    }

    return c_new;
  };

  // Time stepping
  for (int n = 0; n < n_steps; ++n) {
    double t_np1 = (n + 1) * dt;
    double VG_np1 = wf.eval(t_np1);

    // Solve ψ with old charges but new boundary (instantaneous field)
    {
      sol.K = evaluate_kernels(kernels, grid, cfg.species, sol.alpha);
      auto rho = compute_rho_free(grid, cfg.species, sol.K, sol.c, sol.alpha, cfg.redox);
      sol.psi = solve_poisson_dirichlet(grid, sol.K.epsr, rho, 0.0, VG_np1);
    }

    std::vector<double> alpha_new = sol.alpha;
    std::vector<std::vector<double>> c_new = sol.c;

    if (cfg.redox) {
      if (cfg.reaction_first) {
        alpha_new = do_reaction(sol.alpha, sol.c, sol.psi);
      }
    } else {
      alpha_new.assign(nz, 0.0);
    }

    // Kernels for transport (use alpha_new if reaction_first, else current alpha)
    KernelState Kstep = evaluate_kernels(kernels, grid, cfg.species, cfg.reaction_first ? alpha_new : sol.alpha);

    // Transport
    c_new = transport_step(Kstep, sol.psi, sol.alpha, alpha_new, sol.c);

    // Reaction last
    if (cfg.redox && !cfg.reaction_first) {
      alpha_new = do_reaction(sol.alpha, c_new, sol.psi);
    }

    // Update state to t_{n+1}
    sol.alpha = alpha_new;
    sol.c = c_new;

    // Solve ψ self-consistently at t_{n+1}
    sol.K = evaluate_kernels(kernels, grid, cfg.species, sol.alpha);
    auto rho_np1 = compute_rho_free(grid, cfg.species, sol.K, sol.c, sol.alpha, cfg.redox);
    sol.psi = solve_poisson_dirichlet(grid, sol.K.epsr, rho_np1, 0.0, VG_np1);

    // Record
    auto charge = compute_charge(grid, cfg.species, sol.K, sol.c, sol.psi, sol.alpha, cfg.redox, 0.0, VG_np1);
    double Omega = compute_Omega(cfg.T, grid, cfg.species, sol.K, sol.c, sol.psi, sol.alpha, cfg.redox,
                                 cfg.closure_mode).Omega;

    sol.t.push_back(t_np1);
    sol.VG.push_back(VG_np1);
    sol.Q_gate.push_back(charge.Q_gate);
    sol.Q_vol.push_back(charge.Q_vol);
    sol.alpha_bar.push_back(alpha_bar_sites(grid, sol.K, sol.alpha));
    sol.Omega.push_back(Omega);
  }

  // Optional device mapping
  if (cfg.device) {
    sol.ID = map_to_ID(*cfg.device, sol.VG, sol.Q_gate, sol.alpha_bar);
  }

  return sol;
}

void write_dynamic_outputs(const ChannelConfig& cfg,
                           const DynamicSolution& sol,
                           const std::string& subdir) {
  fs::path outdir(cfg.output_dir);
  if (!subdir.empty()) outdir /= subdir;
  ensure_dir(outdir.string());

  // timeseries
  std::vector<std::vector<double>> cols;
  std::vector<std::string> names;
  names = {"t", "VG", "Q_gate", "Q_vol", "alpha_bar", "Omega"};
  cols = {sol.t, sol.VG, sol.Q_gate, sol.Q_vol, sol.alpha_bar, sol.Omega};
  if (!sol.ID.empty()) {
    names.push_back("ID");
    cols.push_back(sol.ID);
  }
  write_table((outdir / "timeseries.dat").string(), names, cols,
              "t [s], VG [V], Q_* [C/m^2], alpha_bar [-], Omega [J/m^2]");

  // final profiles
  // We'll reuse stationary-style output writer by creating a StationarySolution-like shell.
  StationarySolution fin;
  fin.VG = sol.VG.back();
  fin.psi = sol.psi;
  fin.alpha = sol.alpha;
  fin.c = sol.c;
  fin.K = sol.K;
  fin.charge = compute_charge(Grid1D(sol.psi.size(), cfg.d), cfg.species, sol.K, sol.c, sol.psi, sol.alpha, cfg.redox, 0.0, fin.VG);
  std::string fin_subdir = subdir.empty() ? std::string("final") : (subdir + std::string("/final"));
  write_stationary_outputs(cfg, Grid1D(sol.psi.size(), cfg.d), cfg.species, fin, fin_subdir);

  // results.json for dynamic run
  ResultsIndex idx;
  idx.mode = "dynamic";
  idx.config_used = "config_used.ini";
  idx.summary["n_steps"] = std::to_string(sol.t.size() - 1);
  idx.summary["dt"] = std::to_string(cfg.dt);
  idx.summary["t_end"] = std::to_string(cfg.t_end);
  idx.summary["Q_gate_final_C_per_m2"] = std::to_string(sol.Q_gate.back());
  idx.summary["alpha_bar_final"] = std::to_string(sol.alpha_bar.back());

  idx.datasets["timeseries"] = DatasetMeta{"timeseries.dat", names, "Time series outputs"};
  idx.datasets["final_results"] = DatasetMeta{"final/results.json", {"(see final/results.json)"}, "Final profile datasets"};

  write_results_json(outdir.string(), idx);
}

} // namespace channel
