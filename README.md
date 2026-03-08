# CHANNEL

**CHArge aNd ioN NanoscaLe-to-device Link**

CHANNEL is a 1D continuum simulator for charge, ion, and redox physics in nanoscale-to-device workflows. It is designed to consume **PILOTS-generated kernels** and solve the continuum model behind the CHANNEL framework, including electrostatics, ion thermodynamics/transport, optional redox coupling, and an optional OECT-style observable.

## What CHANNEL solves

CHANNEL combines several pieces of physics in one workflow:

- **Poisson electrostatics** with spatially varying dielectric response
- **Multi-species ion equilibrium / transport** in a Poisson–Nernst–Planck-style formulation
- **Scharfetter–Gummel discretization** for drift–diffusion fluxes
- **Optional redox occupancy** with Langmuir / kinetic terms and electrostatic feedback
- **Optional device mapping** to a simple drain-current observable for OECT-style analysis

CHANNEL supports three closure modes:

- **Mode A — Ω-based analytic closure**  
  Uses the standard OMIEC grand-potential closure with kernels such as `ε_r(z;α)`, `h_i(z;α)`, and `Δμ_i^0(z;α)`.
- **Mode B — Ω-based closure with `ω_extra`**  
  Same as Mode A, but adds an optional excess / surrogate energy density `ω_extra(z;α)`.
- **Mode C — μ-closure via `ϕ_ex`**  
  Skips explicit `Ω` and closes the species chemical potentials with per-species excess potentials `ϕ_ex,i(z;α)`.

> CHANNEL does **not** generate nanoscale kernels by itself. Kernel generation is expected to come from **PILOTS** or another upstream workflow.

---

## Build requirements

You need:

- **CMake >= 3.16**
- A **C++17** compiler
- **pkg-config**
- **json-c** development files
- Optional: **OpenMP** for parallel execution

On Linux, install the system dependencies first. Typical package names are similar to:

```bash
# Debian / Ubuntu style example
sudo apt-get install cmake g++ pkg-config libjson-c-dev
```

---

## Build

```bash
git clone <your-fork-or-repo-url>
cd CHANNEL
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

This produces:

- `build/channel` — the CLI executable
- `build/libchannel_core.a` (or platform equivalent) — the internal C++ library target

### Current packaging status

At the moment, the CMake install rule installs the **CLI executable**:

```bash
cmake --install build --prefix ./install
```

If you want CHANNEL to behave like a reusable C++ package, you will likely want to add:

- header installation
- library installation
- exported CMake targets
- package config files

---

## Quick start (CLI)

Run the bundled example:

```bash
./build/channel --config examples/constant_simple.ini
```

Each run writes:

- `results.json` — machine-readable run summary and dataset index
- `*.dat` — plain text columnar outputs

If dynamic mode is enabled, a `dynamic/` directory is also written.

---

## Quick start (Python)

The Python wrapper is intentionally lightweight: it does **not** reimplement the solver. It discovers and launches the compiled `channel` executable, then parses the generated output files.

Install it in editable mode:

```bash
pip install -e python
```

Then use it like this:

```python
import channelio as ch

out = ch.run("examples/constant_simple.ini")
print(out.output_dir)
print(out.stationary)
```

If the executable is not on `PATH`, either:

- build the repo so that `build/channel` exists, or
- set `CHANNEL_EXE=/path/to/channel`, or
- pass `exe="/path/to/channel"` explicitly

### PILOTS-style convenience wrapper

```python
import channelio as ch

out = ch.run_from_pilots(
    pilots_dir="PILOTS_OUT",
    output_dir="CHANNEL_OUT",
    species=[
        dict(name="Na", kappa=+1, r_born=0.18e-9, D=1.33e-9, c_res=1e26),
        dict(name="Cl", kappa=-1, r_born=0.18e-9, D=2.03e-9, c_res=1e26),
    ],
    T=298.15,
    d=100e-9,
    n_cells=200,
    VG=0.2,
    parameterization="A",
    closure_mode="A",
)
```

Useful helper:

```python
kernels = ch.discover_pilots_kernels("PILOTS_OUT", ["Na", "Cl"])
print(kernels)
```

---

## Configuration model

CHANNEL uses an INI configuration file.

### Core sections

- `[general]` — geometry, temperature, run mode, output directory, OpenMP threads
- `[kernel]` — constant kernels or file-based kernel paths
- `[closure]` — closure mode `A`, `B`, or `C`
- `[species.<name>]` — species properties such as valence, diffusion, and reservoir concentration
- `[redox]` — optional redox / counterion coupling
- `[stationary]` — steady-state solve controls
- `[dynamic]` — transient controls and waveform definition
- `[device]` — optional device-level observable mapping
- `[verify]` — thermodynamic / electrostatic consistency checks

### Minimal example

```ini
[general]
T = 298.15
d = 1e-7
n_cells = 200
output_dir = out_example
mode = stationary

[kernel]
source = constant
parameterization = A
epsr_res = 78.0
epsr_const = 78.0
ns_const = 0.0
rho_base_const = 0.0

[species.Na]
kappa = 1
r_born = 1.8e-10
D = 1.33e-9
c_res = 1e23
mobile = true

[species.Cl]
kappa = -1
r_born = 1.8e-10
D = 2.03e-9
c_res = 1e23
mobile = true

[stationary]
enabled = true
VG = 0.01
max_iter = 500
tol = 1e-8
damping = 0.5
enable_feedback = false

[dynamic]
enabled = false

[verify]
enabled = true
```

See `examples/constant_simple.ini` for the bundled version.

---

## Kernel inputs

CHANNEL supports two kernel file formats.

### 1. `Profile1D`

A simple two-column profile on `z`:

```text
# z [m]    value
0.0        78.0
1e-9       78.0
2e-9       77.5
```

### 2. `FieldZAlpha`

A field tabulated on `(z, α)` with an alpha header:

```text
# alpha: 0.0 0.5 1.0
0.0   78.0  60.0  40.0
1e-9  78.0  60.0  40.0
2e-9  77.5  59.0  39.0
```

The solver interpolates these kernels onto the working grid as needed.

### Kernel keys used by CHANNEL

Common kernel inputs include:

- `epsr_file`
- `ns_file`
- `rho_base_file`
- `omega_extra_file` (Mode B)
- `hi_file.<species>`
- `delta_mu0_file.<species>`
- `phi_ex_file.<species>` (Mode C)

---

## Parameterization / double-counting protocol

`[kernel] parameterization` implements the `q0_strategy` convention from the theory notes:

- **A** — explicit `h_i` plus conditional `Δμ_i^0`
- **B** — absorb `h_i` into `Δμ_i^0` and force `h_i ≡ 1`

In the Ω-based formulation, CHANNEL enforces the protocol as:

```text
U_i^ex = U_i^Born + Δμ_i^0
```

with `h_i` appearing only through the ideal / entropic terms.

---

## Outputs

### Stationary outputs

Typical stationary datasets include:

- `psi.dat` — electrostatic potential
- `c_<species>.dat` — concentration profile per species
- `kernels_core.dat` — evaluated core kernels such as `epsr`, `ns`, and `rho_base`
- `results.json` — summary, dataset manifest, and metadata

### Dynamic outputs

When dynamic mode is enabled, CHANNEL additionally writes time-series outputs such as:

- `t`
- `VG`
- `Q_gate`
- `Q_vol`
- `alpha_bar`
- `Omega`
- optional `ID`

---

## Verification and diagnostics

CHANNEL can report internal consistency checks, including:

- Gauss-law / Maxwell charge-closure mismatch
- a Dirichlet-side Maxwell relation check comparing `dΩ/dVG` and `Q_gate`

Enable this through:

```ini
[verify]
enabled = true
```

---

## Notes and limitations

- CHANNEL is currently oriented around **1D continuum problems**.
- Kernel generation is external to this repository.
- The Python wrapper depends on a working `channel` executable.
- The internal C++ library target is built, but the current install step is CLI-first rather than a full exported C++ SDK.

---

## License

This repository is licensed under **GPL-3.0**. See `LICENSE` for details.
