from __future__ import annotations

import configparser
import json
import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple


@dataclass(frozen=True)
class RunOutput:
    """Convenient handle returned by channelio.run()."""
    output_dir: str
    stationary_dir: Optional[str]
    dynamic_dir: Optional[str]
    stationary: Optional[Dict[str, Any]]
    dynamic: Optional[Dict[str, Any]]


def _find_executable(explicit: Optional[str] = None) -> str:
    if explicit:
        p = Path(explicit)
        if not p.exists():
            raise FileNotFoundError(f"CHANNEL executable not found: {explicit}")
        return str(p)

    # 1) environment override
    env = os.environ.get("CHANNEL_EXE")
    if env:
        p = Path(env)
        if p.exists():
            return str(p)

    # 2) local build (repo-relative): build/channel or build/bin/channel
    here = Path(__file__).resolve()
    repo = here.parents[2]  # .../CHANNEL/python/channelio/api.py -> repo root
    candidates = [
        repo / "build" / "channel",
        repo / "build" / "bin" / "channel",
        repo / "channel",
    ]
    for c in candidates:
        if c.exists():
            return str(c)

    # 3) PATH
    from shutil import which
    exe = which("channel")
    if exe:
        return exe

    raise FileNotFoundError(
        "Cannot find CHANNEL executable. Build it with CMake, or set CHANNEL_EXE, or pass exe=..."
    )


def _patch_output_dir(config_path: str, output_dir: str) -> str:
    """
    Create a temporary INI that is identical to config_path but with [general] output_dir overridden.
    Returns the temporary path.
    """
    cp = configparser.ConfigParser()
    cp.read(config_path)
    if "general" not in cp:
        cp["general"] = {}
    cp["general"]["output_dir"] = output_dir

    tmp = tempfile.NamedTemporaryFile("w", suffix=".ini", delete=False)
    with tmp as f:
        cp.write(f)
    return tmp.name


def load_results(output_dir: str) -> Dict[str, Any]:
    """
    Load a results.json (as produced by CHANNEL) from a directory.
    """
    p = Path(output_dir) / "results.json"
    if not p.exists():
        raise FileNotFoundError(f"results.json not found in: {output_dir}")
    return json.loads(p.read_text())


def run(
    config_path: str,
    *,
    output_dir: Optional[str] = None,
    exe: Optional[str] = None,
    extra_args: Optional[Sequence[str]] = None,
    check: bool = True,
) -> RunOutput:
    """
    Run CHANNEL from an INI configuration file.

    Parameters
    ----------
    config_path:
        Path to config.ini
    output_dir:
        If provided, overrides [general] output_dir by writing a temporary config copy.
    exe:
        Path to the CHANNEL executable. If None, auto-discovered.
    extra_args:
        Additional CLI args passed to CHANNEL.
    check:
        If True, raise on non-zero return code.

    Returns
    -------
    RunOutput:
        Includes parsed stationary/dynamic results.json (if present).
    """
    exe_path = _find_executable(exe)

    cfg_to_use = config_path
    tmp_cfg: Optional[str] = None
    if output_dir is not None:
        tmp_cfg = _patch_output_dir(config_path, output_dir)
        cfg_to_use = tmp_cfg

    args = [exe_path, "--config", cfg_to_use]
    if extra_args:
        args.extend(list(extra_args))

    proc = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if check and proc.returncode != 0:
        msg = (
            "CHANNEL run failed.\n"
            f"Command: {' '.join(args)}\n"
            f"Return code: {proc.returncode}\n"
            f"STDOUT:\n{proc.stdout}\n"
            f"STDERR:\n{proc.stderr}\n"
        )
        raise RuntimeError(msg)

    # Derive output directory (either user override, or read from config)
    out_dir = output_dir
    if out_dir is None:
        # best-effort parse
        cp = configparser.ConfigParser()
        cp.read(config_path)
        out_dir = cp.get("general", "output_dir", fallback="out")

    out = Path(out_dir)
    stationary_dir = out / "stationary"
    dynamic_dir = out / "dynamic"

    stationary = load_results(str(stationary_dir)) if stationary_dir.exists() else None
    dynamic = load_results(str(dynamic_dir)) if dynamic_dir.exists() else None

    if tmp_cfg is not None:
        try:
            os.unlink(tmp_cfg)
        except OSError:
            pass

    return RunOutput(
        output_dir=str(out),
        stationary_dir=str(stationary_dir) if stationary_dir.exists() else None,
        dynamic_dir=str(dynamic_dir) if dynamic_dir.exists() else None,
        stationary=stationary,
        dynamic=dynamic,
    )


def _first_existing(base: Path, names: Sequence[str]) -> Optional[str]:
    for n in names:
        p = base / n
        if p.exists():
            return str(p)
    return None


def discover_pilots_kernels(pilots_dir: str, species_names: Sequence[str]) -> Dict[str, str]:
    """
    Best-effort discovery of kernel files in a PILOTS output directory.

    This function is intentionally conservative: it does not assume a specific PILOTS
    schema, but tries common filenames and glob patterns.

    Returns
    -------
    dict with keys:
      - epsr_file, ns_file, rho_base_file
      - hi_file.<sp>, delta_mu0_file.<sp>
    """
    base = Path(pilots_dir)

    out: Dict[str, str] = {}

    out["epsr_file"] = _first_existing(base, [
        "epsr.dat", "epsilon_r.dat", "epsilon_r_profile.dat", "epsr_profile.dat", "epsr.txt"
    ]) or ""

    out["ns_file"] = _first_existing(base, [
        "ns.dat", "n_s.dat", "sites.dat", "ns_profile.dat", "sites_profile.dat"
    ]) or ""

    out["rho_base_file"] = _first_existing(base, [
        "rho_base.dat", "rho_fixed.dat", "rho_base_profile.dat"
    ]) or ""

    # Optional: Mode B ω_extra (if present)
    out["omega_extra_file"] = _first_existing(base, [
        "omega_extra.dat", "omega_extra_profile.dat", "omega_extra.txt"
    ]) or ""

    # Species kernels
    for sp in species_names:
        hi = None
        mu = None

        # common patterns
        for pat in [f"hi_{sp}.dat", f"h_{sp}.dat", f"hi_{sp}.txt", f"h_{sp}.txt"]:
            p = base / pat
            if p.exists():
                hi = str(p)
                break
        if hi is None:
            # glob fallback
            cand = list(base.glob(f"*hi*{sp}*"))
            hi = str(cand[0]) if cand else ""

        for pat in [f"delta_mu0_{sp}.dat", f"mu0_{sp}.dat", f"delta_mu0_{sp}.txt", f"mu0_{sp}.txt"]:
            p = base / pat
            if p.exists():
                mu = str(p)
                break
        if mu is None:
            cand = list(base.glob(f"*mu0*{sp}*")) + list(base.glob(f"*delta_mu0*{sp}*"))
            mu = str(cand[0]) if cand else ""

        out[f"hi_file.{sp}"] = hi or ""
        out[f"delta_mu0_file.{sp}"] = mu or ""

        # Optional: Mode C ϕ_ex
        phi = None
        for pat in [f"phi_ex_{sp}.dat", f"phiex_{sp}.dat", f"phi_ex_{sp}.txt", f"phiex_{sp}.txt"]:
            p = base / pat
            if p.exists():
                phi = str(p)
                break
        if phi is None:
            cand = list(base.glob(f"*phi_ex*{sp}*")) + list(base.glob(f"*phiex*{sp}*"))
            phi = str(cand[0]) if cand else ""
        out[f"phi_ex_file.{sp}"] = phi or ""

    return out


def run_from_pilots(
    pilots_dir: str,
    *,
    output_dir: str,
    species: Sequence[Mapping[str, Any]],
    T: float,
    d: float,
    n_cells: int,
    VG: float = 0.0,
    parameterization: str = "A",
    closure_mode: str = "A",
    epsr_res: float = 78.0,
    stationary: bool = True,
    dynamic: bool = False,
    dynamic_params: Optional[Mapping[str, Any]] = None,
    redox: Optional[Mapping[str, Any]] = None,
    device: Optional[Mapping[str, Any]] = None,
    exe: Optional[str] = None,
) -> RunOutput:
    """
    PILOTS->CHANNEL convenience wrapper.

    It discovers kernel files in `pilots_dir`, writes a CHANNEL config.ini into `output_dir`,
    then runs CHANNEL and returns the parsed results.

    You must still provide `species` meta (valence, D, r_born, c_res) because PILOTS
    outputs are often kernel-only.

    Parameters
    ----------
    pilots_dir:
        Directory produced by PILOTS.
    output_dir:
        CHANNEL output directory (will be created).
    species:
        List of dicts, each with:
          name, kappa, r_born, D, c_res, mobile (optional)
    dynamic_params:
        If dynamic=True, pass dt, t_end, waveform, VG0, VG1, freq, ...
    redox:
        If provided, enables [redox] section.
    device:
        If provided, enables [device] section.
    """
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    sp_names = [str(s["name"]) for s in species]
    kernels = discover_pilots_kernels(pilots_dir, sp_names)

    # Write INI
    ini_path = out / "config.ini"
    cp = configparser.ConfigParser()

    cp["general"] = {
        "T": str(T),
        "d": str(d),
        "n_cells": str(n_cells),
        "output_dir": str(out),
        "mode": "both" if (stationary and dynamic) else ("stationary" if stationary else "dynamic"),
    }

    cp["kernel"] = {
        "source": "files",
        "parameterization": str(parameterization),
        "epsr_res": str(epsr_res),
        "epsr_file": kernels.get("epsr_file", ""),
        "ns_file": kernels.get("ns_file", ""),
        "rho_base_file": kernels.get("rho_base_file", ""),
        "omega_extra_file": kernels.get("omega_extra_file", ""),
    }
    # inline species kernel keys
    for k, v in kernels.items():
        if k.startswith("hi_file.") or k.startswith("delta_mu0_file.") or k.startswith("phi_ex_file."):
            cp["kernel"][k] = v

    cp["closure"] = {"mode": str(closure_mode)}

    cp["stationary"] = {
        "enabled": "true" if stationary else "false",
        "VG": str(VG),
        "max_iter": "2000",
        "tol": "1e-8",
        "damping": "0.4",
        "alpha_init": "1e-6",
        "enable_feedback": "false",
        "explicit_counterion_coupling": "true",
        "compute_capacitance": "false",
    }

    if dynamic:
        dp = dict(dynamic_params or {})
        cp["dynamic"] = {
            "enabled": "true",
            "dt": str(dp.get("dt", 1e-6)),
            "t_end": str(dp.get("t_end", 1e-3)),
            "waveform": str(dp.get("waveform", "step")),
            "VG0": str(dp.get("VG0", VG)),
            "VG1": str(dp.get("VG1", VG)),
            "freq": str(dp.get("freq", 1000.0)),
            "reaction_first": str(dp.get("reaction_first", True)).lower(),
        }
        if "waveform_file" in dp:
            cp["dynamic"]["waveform_file"] = str(dp["waveform_file"])
    else:
        cp["dynamic"] = {"enabled": "false"}

    # Species sections
    for sp in species:
        name = str(sp["name"])
        sec = f"species.{name}"
        cp[sec] = {
            "kappa": str(sp.get("kappa", 0)),
            "r_born": str(sp.get("r_born", 0.0)),
            "D": str(sp.get("D", 0.0)),
            "c_res": str(sp.get("c_res", 0.0)),
            "mobile": str(sp.get("mobile", True)).lower(),
        }

    if redox is not None:
        cp["redox"] = {"enabled": "true"}
        for k, v in redox.items():
            cp["redox"][str(k)] = str(v)
    else:
        cp["redox"] = {"enabled": "false"}

    if device is not None:
        cp["device"] = {"enabled": "true"}
        for k, v in device.items():
            cp["device"][str(k)] = str(v)
    else:
        cp["device"] = {"enabled": "false"}

    cp["verify"] = {"enabled": "true"}

    with open(ini_path, "w") as f:
        cp.write(f)

    return run(str(ini_path), exe=exe, output_dir=str(out))
