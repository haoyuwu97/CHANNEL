#include <channel/profile.hpp>

#include <algorithm>
#include <stdexcept>

namespace channel {

namespace {

void require_monotonic(const std::vector<double>& x, const std::string& name) {
  if (x.size() < 2) {
    throw std::runtime_error(name + ": need at least 2 points");
  }
  for (std::size_t i = 1; i < x.size(); ++i) {
    if (!(x[i] > x[i - 1])) {
      throw std::runtime_error(name + ": grid must be strictly increasing");
    }
  }
}

std::size_t upper_index(const std::vector<double>& x, double xq) {
  // return i such that x[i-1] <= xq < x[i], with i in [1, n-1].
  auto it = std::upper_bound(x.begin(), x.end(), xq);
  if (it == x.begin()) return 1;
  if (it == x.end()) return x.size() - 1;
  return static_cast<std::size_t>(it - x.begin());
}

} // namespace

Profile1D::Profile1D(std::vector<double> z_in, std::vector<double> v_in)
    : z(std::move(z_in)), val(std::move(v_in)) {
  if (z.size() != val.size()) {
    throw std::runtime_error("Profile1D: z and val size mismatch");
  }
  require_monotonic(z, "Profile1D::z");
}

double Profile1D::eval(double z_query) const {
  if (z.empty()) {
    throw std::runtime_error("Profile1D::eval: empty profile");
  }
  if (z_query <= z.front()) return val.front();
  if (z_query >= z.back()) return val.back();

  std::size_t i = upper_index(z, z_query);
  std::size_t i0 = i - 1;
  std::size_t i1 = i;

  double x0 = z[i0], x1 = z[i1];
  double t = (z_query - x0) / (x1 - x0);
  return (1.0 - t) * val[i0] + t * val[i1];
}

std::vector<double> Profile1D::eval_on(const std::vector<double>& z_query) const {
  std::vector<double> out;
  out.reserve(z_query.size());
  for (double zq : z_query) out.push_back(eval(zq));
  return out;
}

std::vector<double> Profile1D::derivative_on_grid() const {
  if (z.size() < 2) return {};

  std::vector<double> d(z.size(), 0.0);
  d.front() = (val[1] - val[0]) / (z[1] - z[0]);
  d.back()  = (val.back() - val[val.size() - 2]) / (z.back() - z[z.size() - 2]);

  for (std::size_t k = 1; k + 1 < z.size(); ++k) {
    d[k] = (val[k + 1] - val[k - 1]) / (z[k + 1] - z[k - 1]);
  }
  return d;
}

// ================= FieldZAlpha =================

FieldZAlpha::FieldZAlpha(std::vector<double> z_in, std::vector<double> alpha_in, std::vector<double> data_in)
    : z(std::move(z_in)), alpha(std::move(alpha_in)), data(std::move(data_in)) {
  require_monotonic(z, "FieldZAlpha::z");
  require_monotonic(alpha, "FieldZAlpha::alpha");
  if (data.size() != z.size() * alpha.size()) {
    throw std::runtime_error("FieldZAlpha: data size mismatch (expected Na*Nz)");
  }
}

std::vector<double> FieldZAlpha::eval_alpha_slice(double alpha_query) const {
  std::vector<double> f;
  std::vector<double> df;
  eval_alpha_slice_with_deriv(alpha_query, f, df);
  return f;
}

void FieldZAlpha::eval_alpha_slice_with_deriv(double alpha_query,
                                              std::vector<double>& f,
                                              std::vector<double>& df_dalpha) const {
  if (z.empty() || alpha.empty()) {
    throw std::runtime_error("FieldZAlpha: empty grid");
  }

  f.assign(z.size(), 0.0);
  df_dalpha.assign(z.size(), 0.0);

  if (alpha.size() == 1) {
    for (std::size_t k = 0; k < z.size(); ++k) {
      f[k] = data[idx(0, k)];
      df_dalpha[k] = 0.0;
    }
    return;
  }

  // clamp
  if (alpha_query <= alpha.front()) alpha_query = alpha.front();
  if (alpha_query >= alpha.back())  alpha_query = alpha.back();

  std::size_t ia = upper_index(alpha, alpha_query);
  std::size_t i0 = ia - 1;
  std::size_t i1 = ia;

  double a0 = alpha[i0], a1 = alpha[i1];
  double t = (alpha_query - a0) / (a1 - a0);
  double inv_da = 1.0 / (a1 - a0);

  for (std::size_t k = 0; k < z.size(); ++k) {
    double f0 = data[idx(i0, k)];
    double f1 = data[idx(i1, k)];
    f[k] = (1.0 - t) * f0 + t * f1;
    df_dalpha[k] = (f1 - f0) * inv_da;
  }
}

} // namespace channel
