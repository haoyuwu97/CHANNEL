#include <channel/grid.hpp>

#include <stdexcept>

namespace channel {

Grid1D::Grid1D(std::size_t n_cells, double length)
    : n_cells_(n_cells), length_(length) {
  if (n_cells_ < 2) {
    throw std::runtime_error("Grid1D: need at least 2 cells");
  }
  if (!(length_ > 0.0)) {
    throw std::runtime_error("Grid1D: length must be positive");
  }
  dz_ = length_ / static_cast<double>(n_cells_);

  zf_.resize(n_cells_ + 1);
  for (std::size_t f = 0; f < zf_.size(); ++f) {
    zf_[f] = dz_ * static_cast<double>(f);
  }

  zc_.resize(n_cells_);
  for (std::size_t j = 0; j < n_cells_; ++j) {
    zc_[j] = (static_cast<double>(j) + 0.5) * dz_;
  }
}

} // namespace channel
