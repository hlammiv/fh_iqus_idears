#!/usr/bin/env python3
"""Extend W_MEASURED past tau = 2, but ONLY if the clamp probe validates.

    python3 calibration/apply_clamp.py            # report, change nothing
    python3 calibration/apply_clamp.py --write    # patch fhcost/hubbard.py

WHY THIS IS A SCRIPT AND NOT A JUDGEMENT CALL
  alpha = 1.75 comes from clamping W at its tau = 2 value, and every plotted
  point sits in the clamped region -- so the clamp, not the measured
  tau-dependence, sets the exponent. The clamp probe measures W beyond tau = 2 on
  TWO lattice sizes, because a small lattice cannot honestly be asked about long
  times: at tau = 4 the cone spans 2 v tau = 16 sites, wider than either patch.

  The decision rule was fixed BEFORE the data arrived:

    sizes AGREE (spread < AGREE_TOL at every tau)
        the tau-trend past 2 is real, not finite-size contamination.
        Extend W_MEASURED with the measured values and report the corrected
        alpha. This is what the user chose over keeping the clamp as a scenario.

    sizes DISAGREE
        a lattice this small cannot test tau > 2. Change nothing; the clamp
        stands for want of evidence, and that is recorded rather than papered
        over.

  Writing the rule down first is the point. It is a one-command update either
  way, and the script refuses to write when the data does not support it.
"""
from __future__ import annotations
import argparse, json, math, pathlib, re, sys

HERE = pathlib.Path(__file__).parent
ROOT = HERE.parent
CLAMP = HERE / "data" / "clamp_probe.json"
AGREE_TOL = 1.5          # max/min across lattice sizes, at every tau


def load():
    """Merge every clamp_probe*.json: the base run plus any single-size reruns.

    A third lattice size can only make the spread test HARDER to pass, never
    easier, so adding one is not a second bite at a rule that already returned
    DISAGREE. The rule itself is unchanged.
    """
    files = sorted((HERE / "data").glob("clamp_probe*.json"))
    if not files:
        sys.exit(f"no clamp data yet at {CLAMP}")
    base = json.loads(files[0].read_text())
    meta, seen = [], set()
    for f in files:
        for m in json.loads(f.read_text())["meta"]:
            k = (m["patch"], m["tau"])
            if k not in seen:
                seen.add(k)
                meta.append(m)
    base["meta"] = meta
    base["sources"] = [f.name for f in files]
    return base


def verdict(C):
    taus = C["clamp_taus"]
    ns = sorted({m["n"] for m in C["meta"]})
    table, spreads = {}, {}
    for tau in taus:
        got = {m["n"]: m["W_eff"] for m in C["meta"] if m["tau"] == tau}
        if len(got) < 2:
            return None, f"only {len(got)} lattice size(s) at tau = {tau}; "\
                         f"the two-size comparison is the whole test"
        table[tau] = got
        spreads[tau] = max(got.values()) / min(got.values())
    agree = all(s < AGREE_TOL for s in spreads.values())
    return {"taus": taus, "ns": ns, "table": table, "spreads": spreads,
            "agree": agree}, None


def main(write=False):
    C = load()
    v, why = verdict(C)
    if v is None:
        print(f"INCOMPLETE: {why}")
        return 1
    sys.path.insert(0, str(ROOT))
    from fhcost import hubbard as H

    print(f"{'tau':>5} " + "".join(f"{'n=' + str(n):>11}" for n in v["ns"])
          + f"{'clamped':>10} {'spread':>8} {'vs clamp':>9}")
    for tau in v["taus"]:
        got, clamped = v["table"][tau], H.w_measured(tau, 4.0)
        # the LARGEST lattice is the least contaminated; use it for the ratio
        best = got[max(got)]
        print(f"{tau:>5} " + "".join(f"{got.get(n, float('nan')):>11.5f}"
                                     for n in v["ns"])
              + f"{clamped:>10.5f} {v['spreads'][tau]:>7.2f}x {best / clamped:>8.2f}x")

    if not v["agree"]:
        print(f"\nVERDICT: the lattice sizes DISAGREE (worst spread "
              f"{max(v['spreads'].values()):.2f}x > {AGREE_TOL}).")
        print("A patch this small cannot be asked about tau > 2 -- at tau = 4 the")
        print("cone is wider than the lattice. The clamp stands for want of")
        print("evidence. NOTHING IS WRITTEN.")
        return 0

    print(f"\nVERDICT: the lattice sizes AGREE (worst spread "
          f"{max(v['spreads'].values()):.2f}x < {AGREE_TOL}).")
    new = {tau: v["table"][tau][max(v["table"][tau])] for tau in v["taus"]}
    ratios = [new[t] / H.w_measured(t, 4.0) for t in v["taus"]]
    print(f"The tau-trend past 2 is real. Measured W is "
          f"{min(ratios):.2f}-{max(ratios):.2f}x the clamped value, so the clamp "
          f"is {'CONSERVATIVE' if max(ratios) < 1 else 'OPTIMISTIC'}.")
    print(f"Step count moves as sqrt(W): "
          f"{math.sqrt(min(ratios)):.2f}-{math.sqrt(max(ratios)):.2f}x.")
    print("\nNew W_MEASURED[4.0] entries: "
          + ", ".join(f"{t}: {w:.4f}" for t, w in new.items()))

    if not write:
        print("\n(report only; pass --write to patch fhcost/hubbard.py)")
        return 0

    src = ROOT / "fhcost" / "hubbard.py"
    text = src.read_text()
    m = re.search(r"( *4\.0: \{[^}]*\},)", text)
    if not m:
        sys.exit("could not find the W_MEASURED U=4 row")
    old = m.group(1)
    entries = dict(re.findall(r"([\d.]+): ([\d.]+)", old))
    for t, w in new.items():
        entries[str(t)] = f"{w:.4f}"
    body = ", ".join(f"{k}: {entries[k]}"
                     for k in sorted(entries, key=float))
    text = text.replace(old, f"    4.0: {{{body}}},")
    # the domain statement has to move with the data -- BOTH axes. Leaving
    # W_DOMAIN behind would make calibration_status report points as
    # out-of-domain that are now inside it, which is the opposite of the bug it
    # was written to catch.
    text = re.sub(r'"tau": \(0\.25, [\d.]+\)',
                  f'"tau": (0.25, {max(v["taus"])})', text)
    text = re.sub(r'"sites": \(4, \d+\)',
                  f'"sites": (4, {max(v["ns"])})', text)
    src.write_text(text)
    print(f"\npatched {src}: U = 4 row extended to tau = {max(v['taus'])}, "
          f"and W_DOMAIN with it")
    print("now run: python3 make_record.py && python3 -m fhcost.selftest")
    return 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--write", action="store_true")
    sys.exit(main(ap.parse_args().write))
