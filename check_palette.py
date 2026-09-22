"""CVD validation for the figure palette. Computed, not eyeballed.

sRGB -> linear -> Machado et al. (2009) dichromat simulation (severity 1.0)
-> OKLab -> Delta E (x100). Targets from the dataviz procedure:
  adjacent-pair CVD Delta E >= 8   (6-8 legal only with secondary encoding)
  normal-vision floor      >= 15   (hard fail below)
"""
import itertools, math

M = {
 "protan":  [[0.152286, 1.052583, -0.204868],
             [0.114503, 0.786281,  0.099216],
             [-0.003882, -0.048116, 1.051998]],
 "deutan":  [[0.367322, 0.860646, -0.227968],
             [0.280085, 0.672501,  0.047413],
             [-0.011820, 0.042940,  0.968881]],
 "tritan":  [[1.255528, -0.076749, -0.178779],
             [-0.078411, 0.930809,  0.147602],
             [0.004733,  0.691367,  0.303900]],
}

def hex2lin(h):
    h = h.lstrip("#")
    c = [int(h[i:i+2], 16)/255 for i in (0, 2, 4)]
    return [x/12.92 if x <= 0.04045 else ((x+0.055)/1.055)**2.4 for x in c]

def apply(mat, v):
    return [max(0.0, min(1.0, sum(mat[i][j]*v[j] for j in range(3)))) for i in range(3)]

def lin2oklab(c):
    r, g, b = c
    l = 0.4122214708*r + 0.5363325363*g + 0.0514459929*b
    m = 0.2119034982*r + 0.6806995451*g + 0.1073969566*b
    s = 0.0883024619*r + 0.2817188376*g + 0.6299787005*b
    l_, m_, s_ = (x ** (1/3) if x > 0 else -((-x) ** (1/3)) for x in (l, m, s))
    return (0.2104542553*l_ + 0.7936177850*m_ - 0.0040720468*s_,
            1.9779984951*l_ - 2.4285922050*m_ + 0.4505937099*s_,
            0.0259040371*l_ + 0.7827717662*m_ - 0.8086757660*s_)

def dE(a, b, mode=None):
    la, lb = hex2lin(a), hex2lin(b)
    if mode:
        la, lb = apply(M[mode], la), apply(M[mode], lb)
    pa, pb = lin2oklab(la), lin2oklab(lb)
    return 100*math.sqrt(sum((x-y)**2 for x, y in zip(pa, pb)))

NAMES = {"#0072B2": "BLUE  NISQ+PEC", "#56B4E9": "SKY   NISQ unmitig.",
         "#D55E00": "ORANGE STAR",    "#009E73": "GREEN surface FT",
         "#333333": "INK   ideal p=0", "#E69F00": "AMBER Pinnacle"}
PAL = list(NAMES)

fails = warns = 0
print(f"{'pair':<44}{'normal':>8}{'protan':>8}{'deutan':>8}{'tritan':>8}  verdict")
for a, b in itertools.combinations(PAL, 2):
    n = dE(a, b)
    cv = [dE(a, b, m) for m in ("protan", "deutan", "tritan")]
    worst = min(cv)
    if n < 15:
        v, fails = "FAIL normal<15", fails+1
    elif worst < 6:
        v, fails = "FAIL cvd<6", fails+1
    elif worst < 8:
        v, warns = "warn cvd 6-8", warns+1
    else:
        v = "pass"
    print(f"{NAMES[a].split()[0]:>7} / {NAMES[b].split()[0]:<8} {'':<26}"[:44]
          + f"{n:>8.1f}" + "".join(f"{x:>8.1f}" for x in cv) + f"  {v}")
print(f"\n{len(list(itertools.combinations(PAL,2)))} pairs: {fails} fail, {warns} warn")
print("all curves are ALSO direct-labelled along the curve (house style: no boxed")
print("legend), so identity never rests on colour alone even where CVD dE is tight.")
raise SystemExit(1 if fails else 0)
