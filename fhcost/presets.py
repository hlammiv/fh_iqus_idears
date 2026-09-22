"""Named assumption sets, so both teams' models can be run from one codebase.

Every field that differs between this repo and an independent parallel resource model
is a knob on Config; these presets just bundle them.

    from fhcost.presets import PRESETS
    cfg = PRESETS["alt_model"]
"""
from __future__ import annotations
from .budget import Config, DEFAULT

# This repo: slide-1 spec (two-time correlator, 10% relative), light cone crossing
# the lattice, Trotter steps from Campbell's commutator bound, 2q noise only.
THIS_REPO = DEFAULT

# The alternative model's central scenario. Values are its declared assumptions,
# not claims about which set is correct.
ARCH_COMPARISON = DEFAULT.but(
    observable="local",          # local doublon density, single-time -> no ancilla
    eps=0.01,                    # 1% relative on D* = 0.1  (eps_abs = 1e-3)
    s_sig=0.1,                   # D* -- the one input we already agreed on
    tmax_mode="const",
    tmax_const=1.0,              # fixed tau, not sqrt(m)/v_B
    trotter_mode="fixed_density",
    steps_per_tau=4.0,           # r = ceil(4 tau), a fixed step density rather than a
                                 # bound-derived count; not calibrated to the accuracy target
    pec_coeff=4.0,               # Gamma^2 = exp(4 p G)
    noise_channels=20.02 / 10.59,  # their chi decomposition: 2q + 1q + idle + SPAM
    star_law="angle",            # 4 alpha p Theta + 4 f R
    star_alpha=1.5,
    star_floor=1e-5,
    magic_source="cultivation",
    lanes=4.0,
    t_round=240e-9,              # 240 ns syndrome cycle
    lightcone_frac=1.0,          # discussion model evolves the full lattice
    routing_power=1.0,           # their compiled anchor scales as L^3, not L^2
)

# Each differing assumption adopted one at a time, for bisecting the gap.
SWAPS = {
    "observable -> local (no ancilla)":  dict(observable="local"),
    "precision -> 1% relative":          dict(eps=0.01),
    "time -> fixed tau = 1":             dict(tmax_mode="const", tmax_const=1.0),
    "Trotter -> r = 4 tau":              dict(trotter_mode="fixed_density"),
    "PEC -> exp(4pG)":                   dict(pec_coeff=4.0),
    "noise -> +1q +idle +SPAM":          dict(noise_channels=20.02 / 10.59),
    "STAR -> angle law":                 dict(star_law="angle"),
    "magic -> cultivation":              dict(magic_source="cultivation"),
    "delivery -> 4 lanes":               dict(lanes=4.0),
    "round -> 240 ns":                   dict(t_round=240e-9),
    "light cone -> off":                 dict(lightcone_frac=1.0),
    "routing -> G ~ L^3 (compiled)":     dict(routing_power=1.0),
}

PRESETS = {"this_repo": THIS_REPO, "alt_model": ARCH_COMPARISON}


if __name__ == "__main__":
    from . import nisq, ftqc, classical
    def row(c):
        fow = c.but(pl_model="fowler")
        lo, hi = classical.band(c)["band"]
        return (nisq.max_m(1e6, c, "pec"), ftqc.max_m_star(1e6, c),
                ftqc.max_m_surface(1e6, fow), hi)

    print("m at n = 1e6, by preset\n")
    print(f"{'preset':<34}{'NISQ+PEC':>9}{'STAR':>8}{'FT':>8}{'classical':>11}")
    for name, c in PRESETS.items():
        a, b, cc, d = row(c)
        print(f"{name:<34}{a:>9.1f}{b:>8.1f}{cc:>8.1f}{d:>11.0f}")

    print("\n\nbisecting the disagreement: one swap at a time from this repo\n")
    base = row(THIS_REPO)
    print(f"{'swap':<34}{'NISQ+PEC':>9}{'STAR':>8}{'FT':>8}")
    print(f"{'(baseline, this repo)':<34}{base[0]:>9.1f}{base[1]:>8.1f}{base[2]:>8.1f}")
    for label, kw in SWAPS.items():
        a, b, cc, _ = row(THIS_REPO.but(**kw))
        print(f"{label:<34}{a:>9.1f}{b:>8.1f}{cc:>8.1f}")
