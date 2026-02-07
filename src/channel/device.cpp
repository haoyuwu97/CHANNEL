#include <channel/device.hpp>

#include <channel/constants.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace channel {

namespace {

double mu_of_alpha(const DeviceParams& dev, double alpha_bar) {
  return dev.mu0 * std::exp(dev.mu_alpha * alpha_bar);
}

} // namespace

std::vector<double> map_to_ID(const DeviceParams& dev,
                              const std::vector<double>& VG,
                              const std::vector<double>& Q_gate,
                              const std::vector<double>& alpha_bar) {
  if (VG.size() != Q_gate.size() || VG.size() != alpha_bar.size()) {
    throw std::runtime_error("map_to_ID: series length mismatch");
  }
  if (dev.Wch <= 0.0 || dev.L <= 0.0) {
    throw std::runtime_error("map_to_ID: Wch and L must be positive");
  }

  std::vector<double> ID(VG.size(), 0.0);

  // Dynamic capacitance estimate if needed
  std::vector<double> Cdyn(VG.size(), 0.0);
  for (std::size_t n = 1; n < VG.size(); ++n) {
    double dV = VG[n] - VG[n - 1];
    double dQ = Q_gate[n] - Q_gate[n - 1];
    if (std::fabs(dV) < 1e-12) {
      Cdyn[n] = Cdyn[n - 1];
    } else {
      Cdyn[n] = dQ / dV;
    }
  }
  if (!Cdyn.empty()) Cdyn[0] = Cdyn.size() > 1 ? Cdyn[1] : 0.0;

  const bool mapA = (dev.map == "A" || dev.map == "a");

  for (std::size_t n = 0; n < VG.size(); ++n) {
    double mu = mu_of_alpha(dev, alpha_bar[n]);
    double Qeff = std::fabs(Q_gate[n]);
    if (mapA) {
      // Use Q* = Cdyn * VG as an alternative proxy (per the CHANNEL theory mapping).
      Qeff = std::fabs(Cdyn[n] * VG[n]);
    }

    // Conductance G ≈ μ * (|Q|) * (Wch / L), current I = G * VD.
    ID[n] = mu * Qeff * (dev.Wch / dev.L) * dev.VD;
  }

  return ID;
}

} // namespace channel
