#!/usr/bin/env python3
"""Verify that `explorer.html` computes the same model as `fhcost/`.

    python3 check_parity.py            # static + headless, exit 1 on any failure
    python3 check_parity.py --static   # skip the browser

Second-pass review finding #3 asked for numerical parity tests between the Python
and JavaScript implementations, across several configurations, on intermediates
as well as maxima. There is no Node on this machine, so the executable half runs
in a headless browser: the explorer's own `parity_check()` recomputes the probe
grid from `RESULTS.json` on load and writes the verdict into the page, and this
script reads it back out of the DOM.

Two layers, because they fail differently:

  STATIC   the generated constants block really is the current Config, and the
           stamped model id really is the current model. Catches a forgotten
           `make_record.py` even if no browser is available.
  HEADLESS the ported FUNCTIONS agree, to 1e-9 on intermediates and 1e-7 on the
           maxima, over every variant in `record.JS_VARIANTS`. A 1e-7 nudge to
           one line of the port trips 209 of these, so it is not a soft check.
"""
from __future__ import annotations
import json, math, pathlib, re, shutil, subprocess, sys
from dataclasses import fields
from fhcost import record
from fhcost.budget import DEFAULT

HERE = pathlib.Path(__file__).resolve().parent
BROWSERS = ("google-chrome", "google-chrome-stable", "chromium", "chromium-browser")


def js_literal(name: str, text: str) -> str:
    """Extract `const NAME = <literal>;` from the generated block."""
    i = text.index(f"const {name} = ") + len(f"const {name} = ")
    depth, j = 0, i
    while True:
        ch = text[j]
        if ch in "{[":
            depth += 1
        elif ch in "}]":
            depth -= 1
        elif ch == ";" and depth == 0:
            return text[i:j]
        j += 1


def static_check() -> list[str]:
    bad = []
    html = (HERE / "explorer.html").read_text()
    rec = json.loads((HERE / "RESULTS.json").read_text())
    live = record.model_id()
    if rec["model_id"] != live:
        bad.append(f"RESULTS.json is model {rec['model_id']}, code is {live}")
    stamped = re.search(r'const MODEL_ID = "([0-9a-f]+)"', html)
    if not stamped:
        bad.append("explorer.html carries no MODEL_ID")
    elif stamped.group(1) != live:
        bad.append(f"explorer.html is model {stamped.group(1)}, code is {live}")
    # every Config field must be present in CFG0 with the same value
    cfg0 = json.loads(js_literal("CFG0", html)
                      .replace("true", "true").replace("false", "false"))
    for f in fields(DEFAULT):
        want = getattr(DEFAULT, f.name)
        if f.name not in cfg0:
            bad.append(f"CFG0 is missing {f.name}")
            continue
        got = cfg0[f.name]
        if isinstance(want, float) and isinstance(got, (int, float)):
            if want != 0 and abs(got - want) / abs(want) > 1e-12:
                bad.append(f"CFG0.{f.name} = {got} vs Config {want}")
        elif got != want:
            bad.append(f"CFG0.{f.name} = {got!r} vs Config {want!r}")
    extra = set(cfg0) - {f.name for f in fields(DEFAULT)}
    if extra:
        bad.append(f"CFG0 has fields Config does not: {sorted(extra)}")
    return bad


def headless_check(timeout_s: int = 240) -> tuple[list[str], str]:
    exe = next((shutil.which(b) for b in BROWSERS if shutil.which(b)), None)
    if exe is None:
        return [], "SKIPPED (no chromium/chrome on PATH)"
    out = subprocess.run(
        [exe, "--headless=new", "--disable-gpu", "--no-sandbox",
         "--virtual-time-budget=60000", "--dump-dom",
         f"file://{HERE / 'explorer.html'}"],
        capture_output=True, text=True, timeout=timeout_s).stdout
    m = re.search(r'id="parity" class="badge (ok|bad)">([^<]*)', out)
    if not m:
        return ["the explorer did not render a parity verdict"], "no verdict"
    return ([] if m.group(1) == "ok" else [m.group(2)]), m.group(2)


def field_sweep_check(timeout_s: int = 600) -> tuple[list[str], str]:
    """Parity over EVERY Config field that changes a probe, not fourteen variants.

    Second-pass #3 left the port "verified, not generated", and the verification
    ran over a hand-picked variant list. A field the JavaScript silently ignores
    passes all of them if none happens to touch it. This builds a throwaway copy
    of the page whose PARITY block is one variant per perturbable field, and
    makes the browser recompute the lot. Nothing shipped grows: the sweep exists
    only for the length of this check.
    """
    exe = next((shutil.which(b) for b in BROWSERS if shutil.which(b)), None)
    if exe is None:
        return [], "SKIPPED (no chromium/chrome on PATH)"
    import make_record
    rec = record.build_field_sweep()
    html = (HERE / "explorer.html").read_text()
    lit = js_literal("PARITY", html)
    sweep = make_record.js({"variants": rec["variants"]})
    # the shipped badge names only the first mismatch, which is useless when the
    # question is WHICH FIELDS the port ignores. Patch it for the throwaway copy.
    badge_old = ('el.textContent=`PARITY MISMATCH (${bad.length}) — believe '
                 'RESULTS.json, not these curves. First: ${bad[0]}`;')
    badge_new = ('el.textContent=`PARITY MISMATCH (${bad.length}) fields: '
                 '${[...new Set(bad.map(x=>x.split("/")[0]))].join(" ")}`;')
    tmp = HERE / ".parity_fields.html"
    page = html.replace(f"const PARITY = {lit};", f"const PARITY = {sweep};", 1)
    if badge_old in page:
        page = page.replace(badge_old, badge_new, 1)
    tmp.write_text(page)
    try:
        out = subprocess.run(
            [exe, "--headless=new", "--disable-gpu", "--no-sandbox",
             "--virtual-time-budget=120000", "--dump-dom", f"file://{tmp}"],
            capture_output=True, text=True, timeout=timeout_s).stdout
    finally:
        tmp.unlink(missing_ok=True)
    m = re.search(r'id="parity" class="badge (ok|bad)">([^<]*)', out)
    n = len(rec["variants"])
    if not m:
        return [f"the explorer did not render a verdict over {n} field variants"], "no verdict"
    msg = (f"{m.group(2)} over {n} perturbable fields; {len(rec['inert'])} more "
           f"move no probe ({', '.join(sorted(rec['inert'])[:6])}...)")
    return ([] if m.group(1) == "ok" else [msg]), msg


def main(argv: list[str]) -> int:
    bad = static_check()
    for b in bad:
        print("[STATIC FAIL]", b)
    if not bad:
        print(f"[static ok] explorer and RESULTS.json are model {record.model_id()}")
    if "--static" not in argv:
        hb, msg = headless_check()
        print(f"[headless] {msg}")
        bad += hb
        if "--fields" in argv:
            fb, fmsg = field_sweep_check()
            print(f"[fields]   {fmsg}")
            bad += fb
    print("PARITY OK" if not bad else f"PARITY FAILED ({len(bad)})")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
