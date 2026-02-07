#pragma once

#include <cstddef>
#include <vector>

namespace channel {

class Grid1D {
public:
  Grid1D() = default;
  Grid1D(std::size_t n_cells, double length);

  std::size_t n_cells() const { return n_cells_; }
  double length() const { return length_; }
  double dz() const { return dz_; }

  // Cell centers z_j (size n_cells)
  const std::vector<double>& zc() const { return zc_; }

  // Face positions z_{j+1/2} (size n_cells+1), including 0 and length.
  const std::vector<double>& zf() const { return zf_; }

private:
  std::size_t n_cells_ = 0;
  double length_ = 0.0;
  double dz_ = 0.0;
  std::vector<double> zc_;
  std::vector<double> zf_;
};

} // namespace channel
