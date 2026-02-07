#include <channel/poisson.hpp>
#include <channel/constants.hpp>

#include <cmath>
#include <stdexcept>

namespace channel {

static void thomas_solve(std::vector<double>& a, // subdiag (a[0] unused)
                         std::vector<double>& b, // diag
                         std::vector<double>& c, // superdiag (c[n-1] unused)
                         std::vector<double>& d  // rhs, becomes solution
) {
  const std::size_t n = b.size();
  for (std::size_t i = 1; i < n; ++i) {
    double m = a[i] / b[i - 1];
    b[i] -= m * c[i - 1];
    d[i] -= m * d[i - 1];
  }
  d[n - 1] /= b[n - 1];
  for (std::size_t i = n - 1; i-- > 0;) {
    d[i] = (d[i] - c[i] * d[i + 1]) / b[i];
  }
}

std::vector<double> solve_poisson_dirichlet(const Grid1D& grid,
                                            const std::vector<double>& epsr,
                                            const std::vector<double>& rho,
                                            double psi0,
                                            double psiD) {
  const std::size_t n = grid.n_cells();
  if (epsr.size() != n || rho.size() != n) {
    throw std::runtime_error("solve_poisson_dirichlet: epsr/rho length mismatch");
  }
  const double dz = grid.dz();
  const double inv_dz2 = 1.0 / (dz * dz);

  std::vector<double> a(n, 0.0);
  std::vector<double> b(n, 0.0);
  std::vector<double> c(n, 0.0);
  std::vector<double> d(n, 0.0);

  for (std::size_t j = 0; j < n; ++j) {
    if (epsr[j] <= 0.0) throw std::runtime_error("solve_poisson_dirichlet: epsr <= 0");

    double epsL = (j == 0) ? epsr[j] : 0.5 * (epsr[j - 1] + epsr[j]);
    double epsR = (j + 1 == n) ? epsr[j] : 0.5 * (epsr[j] + epsr[j + 1]);

    a[j] = -epsilon0 * epsL * inv_dz2;
    b[j] =  epsilon0 * (epsL + epsR) * inv_dz2;
    c[j] = -epsilon0 * epsR * inv_dz2;
    d[j] = rho[j];
  }

  // Apply Dirichlet BCs via ghost-cell elimination
  d[0] -= a[0] * psi0;
  a[0] = 0.0;
  d[n - 1] -= c[n - 1] * psiD;
  c[n - 1] = 0.0;

  // Solve
  thomas_solve(a, b, c, d);
  return d; // solution
}

} // namespace channel
