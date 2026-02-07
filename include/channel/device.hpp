#pragma once

#include <channel/types.hpp>

#include <string>
#include <vector>

namespace channel {

// Simple OECT mapping utilities.
// Given time series of VG(t), Q_gate(t), and alpha_bar(t), compute ID(t).
std::vector<double> map_to_ID(const DeviceParams& dev,
                              const std::vector<double>& VG,
                              const std::vector<double>& Q_gate,
                              const std::vector<double>& alpha_bar);

} // namespace channel
