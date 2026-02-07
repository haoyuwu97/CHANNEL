#include <channel/config.hpp>
#include <channel/constants.hpp>

#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace channel {

namespace {

std::string trim(const std::string& s) {
  std::size_t i = 0;
  while (i < s.size() && std::isspace(static_cast<unsigned char>(s[i]))) ++i;
  std::size_t j = s.size();
  while (j > i && std::isspace(static_cast<unsigned char>(s[j - 1]))) --j;
  return s.substr(i, j - i);
}

std::string to_lower(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return s;
}

bool parse_bool(const std::string& v) {
  std::string x = to_lower(trim(v));
  return (x == "1" || x == "true" || x == "yes" || x == "on");
}

double parse_double(const std::string& v) {
  std::string x = trim(v);
  if (x.empty()) throw std::runtime_error("parse_double: empty");
  char* end = nullptr;
  double out = std::strtod(x.c_str(), &end);
  if (end == x.c_str() || *end != '\0') {
    throw std::runtime_error("parse_double: invalid number '" + v + "'");
  }
  return out;
}

std::size_t parse_size(const std::string& v) {
  double d = parse_double(v);
  if (d < 0.0) throw std::runtime_error("parse_size: negative");
  return static_cast<std::size_t>(d);
}

int parse_int(const std::string& v) {
  double d = parse_double(v);
  return static_cast<int>(d);
}

std::string get_str(const IniMap& ini, const std::string& sec, const std::string& key,
                    const std::string& def) {
  auto sit = ini.find(sec);
  if (sit == ini.end()) return def;
  auto kit = sit->second.find(key);
  if (kit == sit->second.end()) return def;
  return trim(kit->second);
}

std::optional<std::string> get_str_opt(const IniMap& ini, const std::string& sec, const std::string& key) {
  auto sit = ini.find(sec);
  if (sit == ini.end()) return std::nullopt;
  auto kit = sit->second.find(key);
  if (kit == sit->second.end()) return std::nullopt;
  return trim(kit->second);
}

double get_double(const IniMap& ini, const std::string& sec, const std::string& key, double def) {
  auto v = get_str_opt(ini, sec, key);
  return v ? parse_double(*v) : def;
}

int get_int(const IniMap& ini, const std::string& sec, const std::string& key, int def) {
  auto v = get_str_opt(ini, sec, key);
  return v ? parse_int(*v) : def;
}

std::size_t get_size(const IniMap& ini, const std::string& sec, const std::string& key, std::size_t def) {
  auto v = get_str_opt(ini, sec, key);
  return v ? parse_size(*v) : def;
}

bool get_bool(const IniMap& ini, const std::string& sec, const std::string& key, bool def) {
  auto v = get_str_opt(ini, sec, key);
  return v ? parse_bool(*v) : def;
}

bool starts_with(const std::string& s, const std::string& pfx) {
  return s.rfind(pfx, 0) == 0;
}

} // namespace

IniMap parse_ini_file(const std::string& path) {
  std::ifstream f(path);
  if (!f) {
    throw std::runtime_error("Cannot open config INI: " + path);
  }

  IniMap ini;
  std::string current = "general"; // default if no section
  ini[current] = IniSection{};

  std::string line;
  std::size_t lineno = 0;
  while (std::getline(f, line)) {
    ++lineno;

    // Strip comments (# or ;)
    auto hash = line.find('#');
    auto semi = line.find(';');
    std::size_t cut = std::string::npos;
    if (hash != std::string::npos) cut = hash;
    if (semi != std::string::npos) cut = (cut == std::string::npos) ? semi : std::min(cut, semi);
    if (cut != std::string::npos) line = line.substr(0, cut);

    line = trim(line);
    if (line.empty()) continue;

    if (line.front() == '[' && line.back() == ']') {
      current = trim(line.substr(1, line.size() - 2));
      if (current.empty()) {
        throw std::runtime_error("INI parse error: empty section at line " + std::to_string(lineno));
      }
      ini[current] = IniSection{};
      continue;
    }

    auto eq = line.find('=');
    if (eq == std::string::npos) {
      throw std::runtime_error("INI parse error: expected key=value at line " + std::to_string(lineno));
    }
    std::string key = trim(line.substr(0, eq));
    std::string val = trim(line.substr(eq + 1));
    if (key.empty()) {
      throw std::runtime_error("INI parse error: empty key at line " + std::to_string(lineno));
    }
    ini[current][key] = val;
  }

  return ini;
}

ChannelConfig load_config(const std::string& ini_path) {
  ChannelConfig cfg;
  cfg.ini_raw = parse_ini_file(ini_path);

  const IniMap& ini = cfg.ini_raw;

  cfg.T = get_double(ini, "general", "T", cfg.T);
  cfg.d = get_double(ini, "general", "d", cfg.d);
  cfg.n_cells = get_size(ini, "general", "n_cells", cfg.n_cells);
  cfg.output_dir = get_str(ini, "general", "output_dir", cfg.output_dir);
  cfg.mode = to_lower(get_str(ini, "general", "mode", cfg.mode));
  cfg.omp_threads = get_int(ini, "general", "omp_threads", cfg.omp_threads);

  // kernel
  cfg.kernel_source = to_lower(get_str(ini, "kernel", "source", cfg.kernel_source));
  cfg.parameterization = to_lower(get_str(ini, "kernel", "parameterization", cfg.parameterization));
  cfg.epsr_res = get_double(ini, "kernel", "epsr_res", cfg.epsr_res);
  cfg.epsr_const = get_double(ini, "kernel", "epsr_const", cfg.epsr_const);
  cfg.ns_const = get_double(ini, "kernel", "ns_const", cfg.ns_const);
  cfg.rho_base_const = get_double(ini, "kernel", "rho_base_const", cfg.rho_base_const);

  cfg.kernel_files.epsr = get_str(ini, "kernel", "epsr_file", "");
  cfg.kernel_files.ns = get_str(ini, "kernel", "ns_file", "");
  cfg.kernel_files.rho_base = get_str(ini, "kernel", "rho_base_file", "");

  cfg.kernel_files.omega_extra = get_str(ini, "kernel", "omega_extra_file", "");

  // closure
  cfg.closure_mode = to_lower(get_str(ini, "closure", "mode", cfg.closure_mode));

  // stationary
  cfg.stationary_enabled = get_bool(ini, "stationary", "enabled", cfg.stationary_enabled);
  cfg.VG = get_double(ini, "stationary", "VG", cfg.VG);
  cfg.max_iter = get_int(ini, "stationary", "max_iter", cfg.max_iter);
  cfg.tol = get_double(ini, "stationary", "tol", cfg.tol);
  cfg.damping = get_double(ini, "stationary", "damping", cfg.damping);
  cfg.alpha_init = get_double(ini, "stationary", "alpha_init", cfg.alpha_init);
  cfg.enable_feedback = get_bool(ini, "stationary", "enable_feedback", cfg.enable_feedback);
  cfg.explicit_counterion_coupling =
      get_bool(ini, "stationary", "explicit_counterion_coupling", cfg.explicit_counterion_coupling);
  cfg.compute_capacitance = get_bool(ini, "stationary", "compute_capacitance", cfg.compute_capacitance);
  cfg.cap_deltaV = get_double(ini, "stationary", "cap_deltaV", cfg.cap_deltaV);

  // dynamic
  cfg.dynamic_enabled = get_bool(ini, "dynamic", "enabled", cfg.dynamic_enabled);
  cfg.dt = get_double(ini, "dynamic", "dt", cfg.dt);
  cfg.t_end = get_double(ini, "dynamic", "t_end", cfg.t_end);
  cfg.waveform = to_lower(get_str(ini, "dynamic", "waveform", cfg.waveform));
  cfg.VG0 = get_double(ini, "dynamic", "VG0", cfg.VG0);
  cfg.VG1 = get_double(ini, "dynamic", "VG1", cfg.VG1);
  cfg.freq = get_double(ini, "dynamic", "freq", cfg.freq);
  cfg.waveform_file = get_str(ini, "dynamic", "waveform_file", cfg.waveform_file);
  cfg.reaction_first = get_bool(ini, "dynamic", "reaction_first", cfg.reaction_first);

  // verify
  cfg.verify = get_bool(ini, "verify", "enabled", cfg.verify);

  // species
  for (const auto& [sec_name, sec] : ini) {
    if (!starts_with(sec_name, "species.")) continue;
    std::string spname = sec_name.substr(std::string("species.").size());
    if (spname.empty()) continue;

    Species sp;
    sp.name = spname;
    sp.kappa = get_int(ini, sec_name, "kappa", 0);
    sp.q = static_cast<double>(sp.kappa) * eCharge;
    sp.r_born = get_double(ini, sec_name, "r_born", 0.0);
    sp.D = get_double(ini, sec_name, "D", 0.0);
    sp.c_res = get_double(ini, sec_name, "c_res", 0.0);
    sp.mobile = get_bool(ini, sec_name, "mobile", true);
    cfg.species.push_back(sp);

    // kernel species files (if present)
    KernelFileSpec::Sp kspec;
    std::string ksec = "kernel.species." + spname;
    if (ini.find(ksec) != ini.end()) {
      kspec.hi = get_str(ini, ksec, "hi_file", "");
      kspec.delta_mu0 = get_str(ini, ksec, "delta_mu0_file", "");
      kspec.phi_ex = get_str(ini, ksec, "phi_ex_file", "");
      cfg.kernel_files.species[spname] = kspec;
    } else {
      // Also allow inline keys in [kernel]: hi_file.<name>, delta_mu0_file.<name>
      std::string hi_key = "hi_file." + spname;
      std::string mu_key = "delta_mu0_file." + spname;
      std::string ph_key = "phi_ex_file." + spname;
      auto hi_val = get_str_opt(ini, "kernel", hi_key);
      auto mu_val = get_str_opt(ini, "kernel", mu_key);
      auto ph_val = get_str_opt(ini, "kernel", ph_key);
      if (hi_val || mu_val || ph_val) {
        kspec.hi = hi_val ? *hi_val : "";
        kspec.delta_mu0 = mu_val ? *mu_val : "";
        kspec.phi_ex = ph_val ? *ph_val : "";
        cfg.kernel_files.species[spname] = kspec;
      }
    }
  }

  if (cfg.species.empty()) {
    throw std::runtime_error("Config must define at least one [species.<name>] section");
  }

  // redox (optional)
  bool redox_enabled = get_bool(ini, "redox", "enabled", false);
  if (redox_enabled) {
    RedoxParams r;
    r.deltaG0 = get_double(ini, "redox", "deltaG0", 0.0);
    r.KX = get_double(ini, "redox", "KX", 1.0);
    r.sigmaP = get_int(ini, "redox", "sigmaP", 1);
    r.counterion = get_str(ini, "redox", "counterion", "");
    r.koff = get_double(ini, "redox", "koff", 0.0);
    if (r.counterion.empty()) {
      throw std::runtime_error("[redox] enabled, but counterion is empty");
    }
    cfg.redox = r;
  }

  // device (optional)
  bool device_enabled = get_bool(ini, "device", "enabled", false);
  if (device_enabled) {
    DeviceParams d;
    d.map = get_str(ini, "device", "map", d.map);
    d.Wch = get_double(ini, "device", "Wch", d.Wch);
    d.L = get_double(ini, "device", "L", d.L);
    d.d = get_double(ini, "device", "d", cfg.d);
    d.VD = get_double(ini, "device", "VD", d.VD);
    d.mu0 = get_double(ini, "device", "mu0", d.mu0);
    d.mu_alpha = get_double(ini, "device", "mu_alpha", d.mu_alpha);
    cfg.device = d;
  }

  // basic sanity
  if (cfg.T <= 0.0) throw std::runtime_error("T must be positive");
  if (cfg.d <= 0.0) throw std::runtime_error("d must be positive");
  if (cfg.n_cells < 2) throw std::runtime_error("n_cells must be >= 2");
  if (cfg.dt <= 0.0 && cfg.dynamic_enabled) throw std::runtime_error("dynamic dt must be positive");
  if (cfg.t_end <= 0.0 && cfg.dynamic_enabled) throw std::runtime_error("dynamic t_end must be positive");
  if (cfg.damping <= 0.0 || cfg.damping > 1.0) throw std::runtime_error("damping must be in (0,1]");

  // closure sanity
  if (cfg.closure_mode != "a" && cfg.closure_mode != "b" && cfg.closure_mode != "c") {
    throw std::runtime_error("[closure] mode must be one of: A, B, C");
  }

  // q0_strategy sanity (we keep the legacy key name 'parameterization')
  if (cfg.parameterization != "a" && cfg.parameterization != "b") {
    throw std::runtime_error("[kernel] parameterization must be A or B");
  }

  // Mode C requirements: per-species phi_ex must be provided (at least for mobile species).
  if (cfg.closure_mode == "c") {
    if (cfg.enable_feedback) {
      throw std::runtime_error("closure mode C does not support enable_feedback (requires Ω-based closure)");
    }
    if (cfg.explicit_counterion_coupling) {
      throw std::runtime_error("closure mode C does not support explicit_counterion_coupling");
    }
    for (const auto& sp : cfg.species) {
      if (!sp.mobile) continue;
      auto it = cfg.kernel_files.species.find(sp.name);
      std::string ph = (it != cfg.kernel_files.species.end()) ? it->second.phi_ex : "";
      if (ph.empty()) {
        throw std::runtime_error("closure mode C requires phi_ex_file for species '" + sp.name + "'");
      }
    }
  }

  // Mode B: omega_extra is optional; if absent => 0.

  return cfg;
}

} // namespace channel
