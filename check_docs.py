#!/usr/bin/env python3
"""Do the hand-written documents still agree with the model?

    python3 check_docs.py

HEADLINES.md and README.md are GENERATED, so they cannot drift. METHODS.md and
OPEN_ITEMS.md are prose, and prose goes stale silently. In one working session
four typed numbers were found wrong: a 12.9x overcharge that was really 8.8x, a
Trotter step count of 417 that the recalibration had taken to 11.4, an O3
diagnosis that blamed the code family for a magic-engine wall, and a TDVP table
that kept its old reference after the reference changed. Each had survived
because nothing compared it to anything.

This is the comparison. Each entry names a file, a regular expression that
locates the claim and captures the number, and a callable that recomputes it.
A claim that cannot be located is a failure too -- prose gets rewritten, and a
check that quietly stops matching is worse than no check.

Adding an entry is the cost of quoting a model number in prose. That is the
point: it makes the typed number cheap to verify and slightly expensive to write.
"""
from __future__ import annotations
import json, pathlib, re, sys

ROOT = pathlib.Path(__file__).parent
CAL = ROOT / "calibration" / "data"


def _cal(name):
    p = CAL / name
    return json.loads(p.read_text()) if p.exists() else None


def claims():
    """(file, description, regex, expected, rel_tol). regex must capture one number."""
    from fhcost import classical, ftqc, nisq, hubbard, curves
    from fhcost.budget import DEFAULT
    out = []
    PIN = DEFAULT.but(platform=curves.PINNACLE_PLATFORM)
    dec = nisq.experiment_decomposition()

    out += [
        ("METHODS.md", "cone-model overcharge on their circuit",
         r"an \*\*([\d.]+)× overcharge\*\*", dec["overcharge"], 0.02),
        ("METHODS.md", "cone-model Lambda on their circuit",
         r"against \*\*([\d.]+)\*\* from the cone model",
         dec["lambda_model_on_their_circuit"], 0.02),
        ("METHODS.md", "Trotter steps at their operating point",
         r"Trotter step count — ([\d.]+) here vs", dec["steps_model"], 0.02),
        ("METHODS.md", "net Lambda overestimate",
         r"net overestimate of Lambda\*\* \| \*\*([\d.]+)x\*\*",
         dec["net_overestimate"], 0.02),
        ("METHODS.md", "their gates per site per step",
         r"their compiled circuit\n?is ([\d.]+) two-qubit gates per site per step",
         dec["c_g_theirs"], 0.02),
        ("OPEN_ITEMS.md", "cone-model overcharge (O5)",
         r"an \*\*([\d.]+)x overcharge\*\*", dec["overcharge"], 0.02),
        ("OPEN_ITEMS.md", "Pinnacle m ceiling (O3)",
         r"The ceiling is `m ≈ (\d+)`", ftqc.pinnacle_m_ceiling(PIN), 0.02),
        ("OPEN_ITEMS.md", "cleanest published engine p_out (O3)",
         r"`PIN_ENGINE_TABLE` publishes is `1e-(\d+)`",
         -__import__("math").log10(min(r[1] for r in ftqc.PIN_ENGINE_TABLE)), 0.01),
        ("METHODS.md", "free-fermion internal validation",
         r"worst disagreement \*\*([\d.]+) × 10⁻¹⁶\*\*",
         classical.FREE_FERMION["worst_abs_err"] * 1e16, 0.05),
        ("METHODS.md", "free-fermion external validation",
         r"same 7 × 4 instance\n\(`python3 calibration/free_fermion.py --deposit`\): \*\*([\d.]+) × 10⁻¹⁰\*\*",
         classical.FREE_FERMION["deposit_abs_err"] * 1e10, 0.05),
        ("METHODS.md", "pooled TDVP slope",
         r"`err ~ chi\^-([\d.]+)`", abs(classical.TDVP_SLOPE), 0.02),
    ]
    sg = _cal("support_growth.json")
    if sg:
        w = sg["deposit"]["w_mean"]
        out += [
            ("OPEN_ITEMS.md", "mean Heisenberg support over their window",
             r"averages \*\*([\d.]+) of 56\*\*", w, 0.02),
            ("METHODS.md", "mean Heisenberg support over their window",
             r"averages \*\*([\d.]+) of 56 modes\*\*", w, 0.02),
            ("OPEN_ITEMS.md", "measured-support damping fraction",
             r"measured Heisenberg support \| ([\d.]+) \|", w / sg["deposit"]["modes"], 0.02),
            ("OPEN_ITEMS.md", "physical front velocity",
             r"front at weight 10⁻² \| \*\*([\d.]+)\*\*", sg["velocity"]["v_front"][0], 0.02),
            ("OPEN_ITEMS.md", "bulk velocity",
             r"bulk \(mean radius\) \| \*\*([\d.]+)\*\*", sg["velocity"]["v_mean"], 0.02),
            ("OPEN_ITEMS.md", "cluster buffer xi at eps = 1e-2",
             r"\| 10⁻² \| [\d.]+ \| [\d.]+ \| \*\*([\d.]+)\*\*",
             sg["velocity"]["xi_eff"]["1e-02"], 0.03),
            ("OPEN_ITEMS.md", "cluster buffer xi at eps = 1e-12",
             r"\| 10⁻¹² \| [\d.]+ \| [\d.]+ \| \*\*([\d.]+)\*\*",
             sg["velocity"]["xi_eff"]["1e-12"], 0.03),
        ]
    nr = _cal("noise_response.json")
    if nr:
        out += [("OPEN_ITEMS.md", "best-conditioned Lambda from their data",
                 r"attenuation `0\.8096 ± 0\.0039`, \*\*`Lambda = ([\d.]+)`\*\*",
                 nr["lambda_best_point"], 0.02)]
    dc = _cal("domain_check.json")
    if dc and "alpha_measured" in dc:
        out += [
            ("OPEN_ITEMS.md", "trajectory alpha point estimate",
             r"`alpha = ([\d.]+) \+- 0\.\d+`", dc["alpha_measured"], 0.02),
            ("OPEN_ITEMS.md", "trajectory alpha standard error",
             r"`alpha = [\d.]+ \+- ([\d.]+)`", dc["alpha_se"], 0.05),
            ("OPEN_ITEMS.md", "trajectory alpha including n = 4, 5",
             r"flips it to \*\*([\d.]+)\*\*", dc["alpha_all_n"], 0.02),
        ]
    tc = _cal("tdvp_check.json")
    if tc:
        d = tc["dimer"]
        out += [("OPEN_ITEMS.md", "TDVP floor at t = 0.1, chi = 256",
                 r"\| 0\.1 \| 0\.9418 \| ([\d.]+)e-3", d["0.1"]["256"] * 1e3, 0.01)]
    return out


def main():
    fails = []
    n = 0
    for fname, what, rx, expect, tol in claims():
        text = (ROOT / fname).read_text()
        m = re.search(rx, text)
        n += 1
        if not m:
            fails.append(f"  [MISSING] {fname}: {what}\n"
                         f"            pattern {rx!r} no longer matches -- the "
                         f"prose was rewritten and the check went blind")
            continue
        got = float(m.group(1))
        if expect == 0 or abs(got - expect) / abs(expect) > tol:
            fails.append(f"  [STALE]   {fname}: {what}\n"
                         f"            document says {got:g}, model says "
                         f"{expect:g} ({abs(got - expect) / abs(expect):.1%} off)")
        else:
            print(f"  [ok] {fname:14s} {what}  --  {got:g}")
    print()
    if fails:
        print("\n".join(fails))
        print(f"\n{len(fails)} of {n} documented numbers are wrong or unlocatable.")
        return 1
    print(f"all {n} documented numbers match the model")
    return 0


if __name__ == "__main__":
    sys.exit(main())
