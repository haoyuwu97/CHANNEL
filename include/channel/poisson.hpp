#pragma once

#include <channel/grid.hpp>

#include <vector>

namespace channel {

// Solve 1D Poisson with variable dielectric and Dirichlet BCs:
//   -(ε0 ε_r(z) ψ')' = ρ(z),  ψ(0)=psi0, ψ(d)=psiD
// Unknowns ψ are stored at cell centers.
std::vector<double> solve_poisson_dirichlet(const Grid1D& grid,
                                            const std::vector<double>& epsr,    // [nz]
                                            const std::vector<double>& rho,     // [nz]
                                            double psi0,
                                            double psiD);

} // namespace channel
