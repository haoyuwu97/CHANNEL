#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace channel {

// Simple 1D profile on a strictly increasing z grid with linear interpolation.
class Profile1D {
public:
  std::vector<double> z;     // [m]
  std::vector<double> val;   // arbitrary

  Profile1D() = default;
  Profile1D(std::vector<double> z_in, std::vector<double> v_in);

  std::size_t size() const { return z.size(); }

  double eval(double z_query) const;
  std::vector<double> eval_on(const std::vector<double>& z_query) const;

  // dv/dz via centered finite differences (boundaries: one-sided).
  std::vector<double> derivative_on_grid() const;
};

// A z-α field: values tabulated on (alpha_grid, z_grid).
// Minimal: the continuum solver works on a fixed z grid, but α varies;
// so we interpolate along α for each stored z_k (no z interpolation).
class FieldZAlpha {
public:
  std::vector<double> z;        // [m], size Nz
  std::vector<double> alpha;    // [0,1], size Na
  std::vector<double> data;     // row-major: a*Nz + k

  FieldZAlpha() = default;
  FieldZAlpha(std::vector<double> z_in, std::vector<double> alpha_in, std::vector<double> data_in);

  std::size_t nz() const { return z.size(); }
  std::size_t na() const { return alpha.size(); }

  std::vector<double> eval_alpha_slice(double alpha_query) const;

  void eval_alpha_slice_with_deriv(double alpha_query,
                                  std::vector<double>& f,
                                  std::vector<double>& df_dalpha) const;

private:
  std::size_t idx(std::size_t a, std::size_t k) const { return a * nz() + k; }
};

} // namespace channel
