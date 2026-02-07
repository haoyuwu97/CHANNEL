#include <channel/kernels.hpp>
#include <channel/constants.hpp>
#include <channel/io.hpp>

#include <algorithm>
#include <cctype>
#include <fstream>
#include <cmath>
#include <stdexcept>

namespace channel {

namespace {

std::string to_lower(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return s;
}

bool is_field_zalpha_file(const std::string& path) {
  std::ifstream f(path);
  if (!f) return false;
  std::string line;
  for (int k = 0; k < 50 && std::getline(f, line); ++k) {
    if (line.empty()) continue;
    std::string s = line;
    // trim leading spaces
    std::size_t i = 0;
    while (i < s.size() && std::isspace(static_cast<unsigned char>(s[i]))) ++i;
    s = s.substr(i);
    if (s.rfind("#", 0) == 0) s = s.substr(1);
    // trim again
    i = 0;
    while (i < s.size() && std::isspace(static_cast<unsigned char>(s[i]))) ++i;
    s = s.substr(i);
    if (s.rfind("alpha:", 0) == 0) return true;
    if (s.empty()) continue;
    if (s[0] == '#') continue;
  }
  return false;
}

bool z_grids_match(const std::vector<double>& a, const std::vector<double>& b, double tol = 1e-12) {
  if (a.size() != b.size()) return false;
  for (std::size_t i = 0; i < a.size(); ++i) {
    if (std::fabs(a[i] - b[i]) > tol) return false;
  }
  return true;
}

std::size_t upper_index(const std::vector<double>& x, double xq) {
  auto it = std::upper_bound(x.begin(), x.end(), xq);
  if (it == x.begin()) return 1;
  if (it == x.end()) return x.size() - 1;
  return static_cast<std::size_t>(it - x.begin());
}

} // namespace

KernelField KernelField::constant(double value, const Grid1D& grid) {
  Profile1D p(grid.zc(), std::vector<double>(grid.n_cells(), value));
  return KernelField{p, std::nullopt};
}

KernelField KernelField::from_profile(Profile1D p) {
  return KernelField{std::move(p), std::nullopt};
}

KernelField KernelField::from_field(FieldZAlpha f) {
  return KernelField{std::nullopt, std::move(f)};
}

void KernelField::eval(const std::vector<double>& z_query,
                       const std::vector<double>& alpha_field,
                       std::vector<double>& out,
                       std::vector<double>* dout_dalpha) const {
  if (z_query.size() != alpha_field.size()) {
    throw std::runtime_error("KernelField::eval: z_query and alpha_field size mismatch");
  }

  if (profile.has_value()) {
    out = profile->eval_on(z_query);
    if (dout_dalpha) dout_dalpha->assign(z_query.size(), 0.0);
    return;
  }
  if (!field.has_value()) {
    throw std::runtime_error("KernelField::eval: neither profile nor field is set");
  }

  const FieldZAlpha& F = *field;

  // alpha values on F.z grid
  std::vector<double> alpha_on_Fz;
  if (z_grids_match(F.z, z_query)) {
    alpha_on_Fz = alpha_field;
  } else {
    Profile1D a_prof(z_query, alpha_field);
    alpha_on_Fz = a_prof.eval_on(F.z);
  }

  // Evaluate on F.z with local alpha
  std::vector<double> f_Fz(F.z.size(), 0.0);
  std::vector<double> df_Fz(F.z.size(), 0.0);

  if (F.alpha.size() == 1) {
    for (std::size_t k = 0; k < F.z.size(); ++k) {
      f_Fz[k] = F.data[0 * F.nz() + k];
      df_Fz[k] = 0.0;
    }
  } else {
    for (std::size_t k = 0; k < F.z.size(); ++k) {
      double a = alpha_on_Fz[k];
      if (a <= F.alpha.front()) a = F.alpha.front();
      if (a >= F.alpha.back())  a = F.alpha.back();

      std::size_t ia = upper_index(F.alpha, a);
      std::size_t i0 = ia - 1;
      std::size_t i1 = ia;

      double a0 = F.alpha[i0], a1 = F.alpha[i1];
      double t = (a - a0) / (a1 - a0);
      double inv_da = 1.0 / (a1 - a0);

      double f0 = F.data[i0 * F.nz() + k];
      double f1 = F.data[i1 * F.nz() + k];
      f_Fz[k] = (1.0 - t) * f0 + t * f1;
      df_Fz[k] = (f1 - f0) * inv_da;
    }
  }

  if (z_grids_match(F.z, z_query)) {
    out = std::move(f_Fz);
    if (dout_dalpha) *dout_dalpha = std::move(df_Fz);
    return;
  }

  Profile1D f_prof(F.z, f_Fz);
  out = f_prof.eval_on(z_query);
  if (dout_dalpha) {
    Profile1D df_prof(F.z, df_Fz);
    *dout_dalpha = df_prof.eval_on(z_query);
  }
}

std::size_t KernelLibrary::species_index(const std::string& name) const {
  for (std::size_t i = 0; i < sp.size(); ++i) {
    if (sp[i].name == name) return i;
  }
  throw std::runtime_error("KernelLibrary: unknown species '" + name + "'");
}

KernelLibrary build_kernel_library(const ChannelConfig& cfg, const Grid1D& grid) {
  KernelLibrary lib;
  lib.T = cfg.T;
  lib.epsr_res = cfg.epsr_res;
  lib.parameterization = to_lower(cfg.parameterization);
  lib.closure_mode = to_lower(cfg.closure_mode);

  // epsr
  if (cfg.kernel_source == "files") {
    if (!cfg.kernel_files.epsr.empty()) {
      if (is_field_zalpha_file(cfg.kernel_files.epsr)) {
        lib.epsr = KernelField::from_field(read_field_zalpha_file(cfg.kernel_files.epsr));
      } else {
        lib.epsr = KernelField::from_profile(read_profile1d_file(cfg.kernel_files.epsr));
      }
    } else {
      lib.epsr = KernelField::constant(cfg.epsr_const, grid);
    }

    if (!cfg.kernel_files.ns.empty()) {
      lib.ns = KernelField::from_profile(read_profile1d_file(cfg.kernel_files.ns));
    } else {
      lib.ns = KernelField::constant(cfg.ns_const, grid);
    }

    if (!cfg.kernel_files.rho_base.empty()) {
      lib.rho_base = KernelField::from_profile(read_profile1d_file(cfg.kernel_files.rho_base));
    } else {
      lib.rho_base = KernelField::constant(cfg.rho_base_const, grid);
    }

    // Mode B: omega_extra (optional)
    if (!cfg.kernel_files.omega_extra.empty()) {
      if (is_field_zalpha_file(cfg.kernel_files.omega_extra)) {
        lib.omega_extra = KernelField::from_field(read_field_zalpha_file(cfg.kernel_files.omega_extra));
      } else {
        lib.omega_extra = KernelField::from_profile(read_profile1d_file(cfg.kernel_files.omega_extra));
      }
    } else {
      lib.omega_extra = KernelField::constant(0.0, grid);
    }

    // species kernels
    for (const auto& sp : cfg.species) {
      SpeciesKernels sk;
      sk.name = sp.name;
      auto it = cfg.kernel_files.species.find(sp.name);
      std::string hi_file = (it != cfg.kernel_files.species.end()) ? it->second.hi : "";
      std::string mu_file = (it != cfg.kernel_files.species.end()) ? it->second.delta_mu0 : "";

      if (!hi_file.empty()) {
        if (is_field_zalpha_file(hi_file)) sk.hi = KernelField::from_field(read_field_zalpha_file(hi_file));
        else sk.hi = KernelField::from_profile(read_profile1d_file(hi_file));
      } else {
        sk.hi = KernelField::constant(1.0, grid);
      }

      if (!mu_file.empty()) {
        if (is_field_zalpha_file(mu_file)) sk.delta_mu0 = KernelField::from_field(read_field_zalpha_file(mu_file));
        else sk.delta_mu0 = KernelField::from_profile(read_profile1d_file(mu_file));
      } else {
        sk.delta_mu0 = KernelField::constant(0.0, grid);
      }

      // Mode C: phi_ex (optional unless closure_mode==c)
      std::string phi_file = (it != cfg.kernel_files.species.end()) ? it->second.phi_ex : "";
      if (!phi_file.empty()) {
        if (is_field_zalpha_file(phi_file)) sk.phi_ex = KernelField::from_field(read_field_zalpha_file(phi_file));
        else sk.phi_ex = KernelField::from_profile(read_profile1d_file(phi_file));
      } else {
        sk.phi_ex = KernelField::constant(0.0, grid);
      }

      lib.sp.push_back(std::move(sk));
    }
  } else if (cfg.kernel_source == "constant") {
    lib.epsr = KernelField::constant(cfg.epsr_const, grid);
    lib.ns = KernelField::constant(cfg.ns_const, grid);
    lib.rho_base = KernelField::constant(cfg.rho_base_const, grid);
    lib.omega_extra = KernelField::constant(0.0, grid);

    for (const auto& sp : cfg.species) {
      SpeciesKernels sk;
      sk.name = sp.name;
      sk.hi = KernelField::constant(1.0, grid);
      sk.delta_mu0 = KernelField::constant(0.0, grid);
      sk.phi_ex = KernelField::constant(0.0, grid);
      lib.sp.push_back(std::move(sk));
    }
  } else {
    throw std::runtime_error("Unknown kernel.source: " + cfg.kernel_source);
  }

  return lib;
}

KernelState evaluate_kernels(const KernelLibrary& lib,
                             const Grid1D& grid,
                             const std::vector<Species>& species,
                             const std::vector<double>& alpha_field) {
  KernelState st;

  const std::vector<double>& z = grid.zc();
  std::vector<double> depsr;
  lib.epsr.eval(z, alpha_field, st.epsr, &depsr);
  st.depsr_dalpha = std::move(depsr);

  // ns, rho_base are treated as alpha-independent in most workflows,
  // but we still allow them to be KernelFields.
  std::vector<double> tmp;
  lib.ns.eval(z, alpha_field, st.ns, nullptr);
  lib.rho_base.eval(z, alpha_field, st.rho_base, nullptr);

  // Mode B: omega_extra (always available; default 0).
  std::vector<double> dome;
  lib.omega_extra.eval(z, alpha_field, st.omega_extra, &dome);
  st.domega_extra_dalpha = std::move(dome);

  if (lib.sp.size() != species.size()) {
    throw std::runtime_error("evaluate_kernels: kernel species count mismatch with config species");
  }

  st.sp.resize(species.size());

  const bool paramB = (lib.parameterization == "b");

  for (std::size_t i = 0; i < species.size(); ++i) {
    const Species& sp = species[i];
    const SpeciesKernels& sk = lib.sp[i];
    if (sk.name != sp.name) {
      throw std::runtime_error("evaluate_kernels: species order mismatch (kernel vs config)");
    }

    auto& s = st.sp[i];

    std::vector<double> dhi;
    sk.hi.eval(z, alpha_field, s.hi, &dhi);

    std::vector<double> dmu0;
    sk.delta_mu0.eval(z, alpha_field, s.delta_mu0, &dmu0);
    s.ddelta_mu0_dalpha = std::move(dmu0);

    // parameterization B => force h_i(z) = 1, ignore its α-derivative.
    if (paramB) {
      s.hi.assign(z.size(), 1.0);
      s.dlnhi_dalpha.assign(z.size(), 0.0);
    } else {
      s.dlnhi_dalpha.resize(z.size(), 0.0);
      for (std::size_t k = 0; k < z.size(); ++k) {
        double h = s.hi[k];
        if (h <= 0.0) throw std::runtime_error("h_i(z) <= 0 encountered in kernels");
        s.dlnhi_dalpha[k] = dhi[k] / h;
      }
    }

    // Born + short-range
    s.Uex.assign(z.size(), 0.0);
    s.dUex_dalpha.assign(z.size(), 0.0);

    for (std::size_t k = 0; k < z.size(); ++k) {
      double eps = st.epsr[k];
      if (eps <= 0.0) throw std::runtime_error("epsilon_r(z) <= 0 encountered in kernels");

      double born = 0.0;
      double dborn = 0.0;
      if (sp.r_born > 0.0 && sp.q != 0.0) {
        double pref = (sp.q * sp.q) / (8.0 * pi * epsilon0 * sp.r_born);
        born = pref * (1.0 / eps - 1.0 / lib.epsr_res);
        dborn = pref * (-1.0) * st.depsr_dalpha[k] / (eps * eps);
      }

      // q0_strategy A/B: Uex includes Born + Δμ0 only.
      // The accessible-volume factor h_i enters only through the ideal term
      // (ln(c/(h c_res))) and through the drift potential (−ln h).
      s.Uex[k] = born + s.delta_mu0[k];
      s.dUex_dalpha[k] = dborn + s.ddelta_mu0_dalpha[k];
    }

    // Mode C: phi_ex and derivative (dimensionless). Always evaluate (cheap); default 0.
    std::vector<double> dphi;
    sk.phi_ex.eval(z, alpha_field, s.phi_ex, &dphi);
    s.dphi_ex_dalpha = std::move(dphi);
  }

  return st;
}

} // namespace channel
