#pragma once

#include <channel/types.hpp>

#include <cstddef>
#include <map>
#include <optional>
#include <string>
#include <vector>

namespace channel {

// Parsed INI as section->(key->value).
using IniSection = std::map<std::string, std::string>;
using IniMap = std::map<std::string, IniSection>;

struct KernelFileSpec {
  // Common profiles for the continuum.
  std::string epsr;      // Profile1D or FieldZAlpha file
  std::string ns;        // Profile1D file (sites number density)
  std::string rho_base;  // Profile1D file

  // Optional extra energy density term ω_extra(z;α) [J/m^3] for closure Mode B.
  // If empty, treated as 0.
  std::string omega_extra;

  // Map species name -> file spec for its kernels
  struct Sp {
    std::string hi;           // Profile1D or FieldZAlpha
    std::string delta_mu0;    // Profile1D or FieldZAlpha, [J]

    // Optional μ-closure potential ϕ_ex(z;α) [dimensionless] for closure Mode C.
    // If empty and Mode C is selected, CHANNEL will error out.
    std::string phi_ex;
  };
  std::map<std::string, Sp> species;
};

struct ChannelConfig {
  // [general]
  double T = 298.15;                // [K]
  double d = 100e-9;                // [m]
  std::size_t n_cells = 200;
  std::string output_dir = "out";
  std::string mode = "stationary";  // stationary | dynamic | both
  int omp_threads = 0;              // 0 => leave as-is

  // [kernel]
  std::string kernel_source = "constant"; // constant | files
  // q0_strategy: "A" (explicit h, conditional Δμ0) or "B" (lumped Δμ0, force h=1)
  std::string parameterization = "A";
  double epsr_res = 78.0;                 // reservoir dielectric
  double epsr_const = 78.0;               // constant kernel fallback
  double ns_const = 0.0;
  double rho_base_const = 0.0;
  KernelFileSpec kernel_files;            // if kernel_source == files

  // [closure]
  // A: analytic Ω-based closure (default)
  // B: Ω-based with additional ω_extra(z;α)
  // C: μ-closure via per-species ϕ_ex(z;α)
  std::string closure_mode = "A";

  // Species list: each in [species.<name>]
  std::vector<Species> species;

  // [redox] (optional)
  std::optional<RedoxParams> redox;

  // [stationary]
  bool stationary_enabled = true;
  double VG = 0.0;                 // [V] (ψ(0)=0, ψ(d)=VG)
  int max_iter = 2000;
  double tol = 1e-8;
  double damping = 0.4;            // Picard mixing for ψ and α
  double alpha_init = 1e-6;
  bool enable_feedback = false;    // include α-feedback terms in δΩ/δα
  bool explicit_counterion_coupling = true; // include ω_redox dependence on cX
  bool compute_capacitance = false;
  double cap_deltaV = 1e-3;        // [V] finite difference for C*_eq

  // [dynamic]
  bool dynamic_enabled = false;
  double dt = 1e-6;               // [s]
  double t_end = 1e-3;            // [s]
  std::string waveform = "step";  // step | sine | triangle | file
  double VG0 = 0.0;
  double VG1 = 0.1;
  double freq = 1000.0;           // [Hz] for periodic waveforms
  std::string waveform_file = "";
  bool reaction_first = true;     // operator splitting order: reaction then transport

  // [device] (optional)
  std::optional<DeviceParams> device;

  // [verify]
  bool verify = true;

  // Convenience: report what was parsed.
  IniMap ini_raw;
};

IniMap parse_ini_file(const std::string& path);
ChannelConfig load_config(const std::string& ini_path);

} // namespace channel
