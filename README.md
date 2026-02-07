# CHANNEL — CHArge aNd ioN NanoscaLe-to-device Link

A continuum **charge + ion** simulator designed to take **PILOTS**-processed nanoscale kernel outputs
and solve the 1D continuum theory described in the CHANNEL framework:

- electrostatics (Poisson) with spatially varying dielectric
- multi-species ion thermodynamics and transport (Poisson–Nernst–Planck with Scharfetter–Gummel)
- optional redox (Langmuir / kinetics) + electrostatic feedback on redox energetics
- optional mapping to a simple OECT drain current observable

CHANNEL supports three **closure modes** (vNext):

* **Mode A (Ω-based, analytic)**: the default OMIEC grand-potential closure using
  `ε_r(z;α)`, `h_i(z;α)`, `Δμ_i^0(z;α)`, ...
* **Mode B (Ω-based + ω_extra)**: identical to Mode A but with an additional, optional
  energy density kernel `ω_extra(z;α)` for unknown/excess physics (still auditable).
* **Mode C (μ-closure via ϕ_ex)**: skip explicit Ω; provide per-species
  `ϕ_ex,i(z;α)` so that `μ_i/(kBT)=ln(c_i/c_i^{res})+ϕ_ex,i+βq_iψ`.

This repository provides:

- `libchannel` (C++ library)
- `channel` (CLI executable)
- `python/channelio` (Python wrapper with PILOTS-style ergonomics)

> Note: Kernel generation (ε_r(z,α), h_i(z,α), Δμ_i^0(z,α), n_s(z), ρ_base(z), …) is expected from PILOTS.

---

## Build (C++)

```bash
cd CHANNEL
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

The executable will be `build/channel`.

---

## Run (CLI)

```bash
./build/channel --config examples/constant_simple.ini
```

Outputs are written under `[general] output_dir`, into subfolders:

- `output_dir/stationary/*`
- `output_dir/dynamic/*` (if enabled)

Each run writes `results.json` plus plain text `.dat` tables.

---

## Python interface

Install the python wrapper in editable mode:

```bash
pip install -e python
```

Then:

```python
import channelio as ch

out = ch.run("config.ini")
print(out.stationary)  # parsed results.json from output_dir/stationary
```

PILOTS interop helper:

```python
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

If you have a custom PILOTS file schema, use `channelio.discover_pilots_kernels()` and write your own config.

---

## Configuration overview

CHANNEL uses an INI config file. Key sections:

- `[general]`: geometry, temperature, output dir, mode
- `[kernel]`: constant or file-based kernels
- `[closure]`: closure mode (A/B/C)
- `[species.<name>]`: ion properties
- `[redox]`: optional redox + counterion coupling
- `[stationary]`: steady-state solver parameters
- `[dynamic]`: transient solver parameters
- `[device]`: optional OECT mapping parameters
- `[verify]`: enable thermodynamic consistency checks

See `examples/constant_simple.ini` for a minimal config.

### q0_strategy / double-counting protocol

`[kernel] parameterization` is the **q0_strategy** from the theory document:

* `A`: explicit `h_i` (accessible volume) and **conditional** `Δμ_i^0` (beyond sieving)
* `B`: lump `h_i` into `Δμ_i^0` and force `h_i ≡ 1`

CHANNEL enforces this protocol by construction:
`U_i^{ex} = U_i^{Born} + Δμ_i^0` and `h_i` enters only through `ln(c/(h c_res))` and `-ln h`.

### Mode B / Mode C extra kernel files

* `[kernel] omega_extra_file`: optional `Profile1D` or `FieldZAlpha` for `ω_extra(z;α)` [J/m^3]
* `[kernel] phi_ex_file.<name>` or `[kernel.species.<name>] phi_ex_file`: optional `Profile1D` or
  `FieldZAlpha` for `ϕ_ex,i(z;α)` [dimensionless] (required for Mode C)

---

## Kernel file formats

**Profile1D** (`z value` per line, comment lines start with `#` or `;`):

```
# z [m]   epsr [-]
0.0      78.0
1e-9     78.0
...
```

**FieldZAlpha** (α header + data rows; each row is `z f(α0) f(α1) ...`):

```
# alpha: 0.0 0.5 1.0
0.0   78.0  60.0  40.0
1e-9  78.0  60.0  40.0
...
```

---

## License

Internal/research use. Add your preferred license here.
