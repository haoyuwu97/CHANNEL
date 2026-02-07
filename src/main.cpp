#include <channel/config.hpp>
#include <channel/grid.hpp>
#include <channel/io.hpp>
#include <channel/kernels.hpp>
#include <channel/stationary.hpp>
#include <channel/dynamic.hpp>
#include <channel/verify.hpp>

#include <iostream>
#include <filesystem>
#include <stdexcept>
#include <string>

#ifdef CHANNEL_HAS_OPENMP
#include <omp.h>
#endif

namespace {

void print_usage() {
  std::cout << "CHANNEL (CHArge aNd ioN NanoscaLe-to-device Link)\n"
            << "Usage:\n"
            << "  channel --config <path/to/config.ini>\n";
}

} // namespace

int main(int argc, char** argv) {
  try {
    std::string cfg_path;
    for (int i = 1; i < argc; ++i) {
      std::string a = argv[i];
      if (a == "--config" && i + 1 < argc) {
        cfg_path = argv[++i];
      } else if (a == "-h" || a == "--help") {
        print_usage();
        return 0;
      } else {
        std::cerr << "Unknown argument: " << a << "\n";
        print_usage();
        return 2;
      }
    }
    if (cfg_path.empty()) {
      print_usage();
      return 2;
    }

    channel::ChannelConfig cfg = channel::load_config(cfg_path);

#ifdef CHANNEL_HAS_OPENMP
    if (cfg.omp_threads > 0) {
      omp_set_num_threads(cfg.omp_threads);
    }
#endif

    // Create output root, copy config
    channel::ensure_dir(cfg.output_dir);
    channel::copy_file(cfg_path, (std::filesystem::path(cfg.output_dir) / "config_used.ini").string());

    channel::Grid1D grid(cfg.n_cells, cfg.d);
    channel::KernelLibrary kernels = channel::build_kernel_library(cfg, grid);

    bool do_stationary = cfg.stationary_enabled && (cfg.mode == "stationary" || cfg.mode == "both");
    bool do_dynamic = cfg.dynamic_enabled && (cfg.mode == "dynamic" || cfg.mode == "both");

    if (!do_stationary && !do_dynamic) {
      std::cerr << "Nothing to do: check [general] mode and [stationary]/[dynamic] enabled flags.\n";
      return 3;
    }

    if (do_stationary) {
      auto sol = channel::solve_stationary(cfg, grid, kernels, cfg.VG, true);
      channel::write_stationary_outputs(cfg, grid, cfg.species, sol, "stationary");

      if (cfg.verify) {
        auto rep = channel::verify_stationary(cfg, grid, kernels, sol);
        std::cout << "[verify] Maxwell mismatch: " << rep.maxwell_mismatch
                  << " (ok=" << (rep.maxwell_ok ? "true" : "false") << ")\n";
        std::cout << "[verify] dOmega/dVG: " << rep.dOmega_dVG
                  << " vs Q_gate: " << rep.Q_gate
                  << " rel_err=" << rep.energy_rel_error
                  << " (ok=" << (rep.energy_ok ? "true" : "false") << ")\n";
      }
    }

    if (do_dynamic) {
      auto sol = channel::run_dynamic(cfg, grid, kernels);
      channel::write_dynamic_outputs(cfg, sol, "dynamic");
    }

    std::cout << "Done. Output in: " << cfg.output_dir << "\n";
    return 0;
  } catch (const std::exception& e) {
    std::cerr << "ERROR: " << e.what() << "\n";
    return 1;
  }
}
