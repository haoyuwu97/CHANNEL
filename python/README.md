# channelio

Python wrapper for the **CHANNEL** continuum simulator.

Typical usage:

```python
import channelio as ch

out = ch.run("config.ini")
print(out.stationary)  # parsed results.json (if produced)
```

PILOTS interop (best-effort discovery of kernel files):

```python
out = ch.run_from_pilots(
    pilots_dir="PILOTS_OUT",
    output_dir="CHANNEL_OUT",
    species=[
        dict(name="Na", kappa=+1, r_born=0.18e-9, D=1.33e-9, c_res=1000*6.022e23),
        dict(name="Cl", kappa=-1, r_born=0.18e-9, D=2.03e-9, c_res=1000*6.022e23),
    ],
    T=298.15,
    d=100e-9,
    n_cells=200,
    VG=0.2,
    parameterization="A",
)
```
