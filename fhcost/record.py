"""The SHARED RESULT RECORD: one executable model, one serialised answer.

Second-pass review finding #3 was that the Python package, `explorer.html`, and
the prose documents reported different numbers at the same nominal settings --
not to rounding, but because they were different models. Three of them had been
maintained by hand.

This module makes the Python package the single source and everything else a
CONSUMER of a versioned record:

    python3 make_record.py            # writes RESULTS.json, stamps explorer.html
    python3 make_record.py --check    # non-zero exit if anything is stale

**model_id** is a hash of the model source AND the default configuration, so any
edit to either invalidates every stamped artefact. The figure, the crossover
table, the headline file, and the explorer all carry it; if two of them show
different ids they were produced by different models and must not be compared.

The record deliberately carries INTERMEDIATES, not only the maxima. Two models
can agree on `max_m` at one n by cancelling errors; they cannot agree on steps,
signal, gate counts, shot counts, code distance, workspace, factory choice,
footprint and per-shot time at a dozen operating points unless they are the same
model. The probe grid includes zero-capacity points and factory/code transitions
for the same reason.
"""
from __future__ import annotations
import hashlib, json, math, pathlib
from dataclasses import fields, asdict
from .budget import Config, DEFAULT
from . import hubbard, nisq, ftqc, classical, curves

SRC = ("budget.py", "hubbard.py", "nisq.py", "ftqc.py", "classical.py",
       "converged.py", "curves.py", "record.py")

# Probe configurations. Each is a NAMED SCENARIO: the explorer must reproduce
# every one of them, so a knob that exists in only one of the two is caught.
VARIANTS = {
    "default":      {},
    "p_1e-4":       {"p": 1e-4},
    "eps_0.01":     {"eps": 0.01},
    "tmax_const":   {"tmax_mode": "const"},
    "tmax_inv_m":   {"tmax_mode": "inv_m"},
    "mpf_k1":       {"trotter_order_k": 1},
    "mpf_k3":       {"trotter_order_k": 3},
    "trotter_bound": {"trotter": "extensive"},
    "u8":           {"U_over_J": 8.0},
    "u0":           {"U_over_J": 0.0},
    "willow":       {"pl_model": "willow"},
    "sig_fixed":    {"signal_regime": "fixed"},
    "budget_day":   {"budget_s": 86400.0},
    # the independent parallel model's declared assumption set, as a whole. If the
    # explorer's "alt model" preset ever drifts from fhcost.presets, this catches it.
    "arch_comparison": None,        # filled in below from the preset
    # hardware platforms: connectivity, gate parallelism and clock
    "helios": None,
    "neutral_atom": None,
    "two_layer_sc": None,
}
def _preset_overrides(name: str) -> dict:
    from dataclasses import fields as _fields
    from .presets import PRESETS
    c = PRESETS[name]
    return {f.name: getattr(c, f.name) for f in _fields(c)
            if getattr(c, f.name) != getattr(DEFAULT, f.name)}


M_PROBES = (4.0, 9.0, 16.0, 36.0, 64.0, 256.0)
N_PROBES = (1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e10)


def model_id() -> str:
    """sha256 over the model source and the default config. 12 hex chars."""
    h = hashlib.sha256()
    here = pathlib.Path(__file__).resolve().parent
    for name in SRC:
        h.update(name.encode())
        h.update((here / name).read_bytes())
    for f in sorted(fields(DEFAULT), key=lambda x: x.name):
        h.update(f"{f.name}={getattr(DEFAULT, f.name)!r};".encode())
    return h.hexdigest()[:12]


def _f(x) -> float | None:
    """JSON-safe float: inf and nan become null so a consumer cannot mistake them."""
    if x is None:
        return None
    x = float(x)
    return x if math.isfinite(x) else None


def circuit_probe(m: float, cfg: Config) -> dict:
    """Everything between the lattice size and the architecture models."""
    c = hubbard.counts(m, cfg)
    return {
        "t_max":      _f(c["t_max"]),
        "signal":     _f(hubbard.signal_at(m, cfg)),
        "eps_abs":    _f(hubbard.eps_absolute(cfg, m)),
        "steps":      _f(c["steps"]),
        "damp_frac":  _f(c["damp_frac"]),
        "g_total":    _f(c["g_total"]),
        "g_cone":     _f(c["g_cone"]),
        "n_rot":      _f(c["n_rot"]),
        "n_rot_cone": _f(c["n_rot_cone"]),
        "depth":      _f(c["depth"]),
        "t_circuit":  _f(c["t_circuit"]),
        "q_per_copy": _f(c["q_per_copy"]),
        "lambda":     _f(nisq.lambda_of(m, cfg)),
        "log_gamma2": _f(nisq.log_gamma_sq(m, cfg)),
        "n_shots":    _f(ftqc.n_shots_total(cfg, m)),
        "t_nisq_pec": _f(nisq.time_required(m, cfg, "pec")),
    }


def arch_probe(m: float, cfg: Config) -> dict:
    """Per-architecture intermediates -- distance, workspace, factory, footprint."""
    out = {}
    s = ftqc.surface_point(m, cfg)
    out["surface"] = None if s is None else {
        "d": s["d"], "q_L": _f(s["q_L"]), "n_t": _f(s["n_t"]), "d_t": _f(s["d_t"]),
        "rounds": _f(s["rounds"]), "n_fac": _f(s["n_fac"]), "phys": _f(s["phys"]),
        "factory": s["factory"], "workspace": _f(s["workspace"]),
        "t_shot": _f(s["t_shot"])}
    st = ftqc.star_point(m, cfg)
    out["star"] = None if st is None else {
        "d": st["d"], "q_L": _f(st["q_L"]), "lam": _f(st["lam"]),
        "phys": _f(st["phys"]), "t_shot": _f(st["t_shot"])}
    pn = ftqc.pinnacle_point(m, cfg)
    out["pinnacle"] = None if pn is None else {
        "d": pn["d"], "k": pn["k"], "blocks": pn["blocks"],
        "cycles": _f(pn["cycles"]), "phys": _f(pn["phys"]),
        "t_shot": _f(pn["t_shot"])}
    return out


def reach_probe(n: float, cfg: Config) -> dict:
    """The maxima. Zeros are kept: a zero IS the answer at small n."""
    return {
        "ideal":     _f(nisq.max_m_ideal(n, cfg)),
        "nisq_none": _f(nisq.max_m(n, cfg, "none")),
        "nisq_pec":  _f(nisq.max_m(n, cfg, "pec")),
        "nisq_zne2": _f(nisq.max_m(n, cfg, "zne2")),
        "star":      _f(ftqc.max_m_star(n, cfg)),
        "surface":   _f(ftqc.max_m_surface(n, cfg)),
        # Pinnacle is costed on the hardware its codes require, exactly as
        # curves.evaluate plots it -- otherwise the record and the figure would
        # disagree, which is the whole point of having one record.
        "pinnacle":  _f(ftqc.max_m_pinnacle(
            n, cfg.but(platform=curves.PINNACLE_PLATFORM))),
        "pinnacle_on_grid": _f(ftqc.max_m_pinnacle(n, cfg)),
    }


def build() -> dict:
    out = {
        "model_id": model_id(),
        "note": "Generated by make_record.py. Do not hand-edit. Every consumer "
                "(explorer.html, README.md, HEADLINES.md, figures) carries this "
                "model_id; differing ids mean differing models.",
        "config": {f.name: getattr(DEFAULT, f.name) for f in fields(DEFAULT)},
        "m_probes": list(M_PROBES),
        "n_probes": list(N_PROBES),
        "variants": {},
    }
    for name, over in VARIANTS.items():
        if over is None:
            over = _preset_overrides({"arch_comparison": "alt_model"}.get(name, name))
        cfg = DEFAULT.but(**over) if over else DEFAULT
        b = classical.band(cfg)
        out["variants"][name] = {
            "overrides": over,
            "circuit": {f"{m:g}": circuit_probe(m, cfg) for m in M_PROBES},
            "arch":    {f"{m:g}": arch_probe(m, cfg) for m in M_PROBES},
            "reach":   {f"{n:g}": reach_probe(n, cfg) for n in N_PROBES},
            "classical_band": [_f(b["band"][0]), _f(b["band"][1])],
        }
    return out


def write(path: str = "RESULTS.json") -> dict:
    rec = build()
    pathlib.Path(path).write_text(json.dumps(rec, indent=1, sort_keys=False) + "\n")
    return rec


def load(path: str = "RESULTS.json") -> dict:
    return json.loads(pathlib.Path(path).read_text())


def compare(a: dict, b: dict, rtol: float = 1e-6) -> list[str]:
    """Recursive numeric diff. Returns a list of human-readable disagreements."""
    out = []

    def walk(x, y, path):
        if isinstance(x, dict) and isinstance(y, dict):
            for k in sorted(set(x) | set(y)):
                if k not in x or k not in y:
                    out.append(f"{path}.{k}: present in only one record")
                else:
                    walk(x[k], y[k], f"{path}.{k}")
        elif isinstance(x, (int, float)) and isinstance(y, (int, float)) \
                and not isinstance(x, bool) and not isinstance(y, bool):
            den = max(abs(float(x)), abs(float(y)), 1e-300)
            if abs(float(x) - float(y)) / den > rtol:
                out.append(f"{path}: {x!r} vs {y!r}")
        elif x != y:
            out.append(f"{path}: {x!r} vs {y!r}")

    walk(a, b, "")
    return out
