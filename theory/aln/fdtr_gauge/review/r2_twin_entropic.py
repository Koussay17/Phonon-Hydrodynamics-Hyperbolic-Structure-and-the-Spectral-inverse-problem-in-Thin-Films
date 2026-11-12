"""Cross-check: 1D closed form with 'local' / 'entropic' interface temperature against src/twin_fdtr.kernel at k = 0,
for a GK AlN film directly on sapphire (no transducer), and rank of the 4-parameter Jacobian on the resonance in SI units.
"""
import sys, json
import numpy as np
sys.path.insert(0, r"C:\Users\Koussay\these\src")
import twin_fdtr as tw

C_aln = 2.41e6
lam = 188.6
e = 500e-9
sap = tw.Material(35.0, 3.03e6)


def closed(freq, lam, C, e, tR, tl, conv, R=0.0):
    p = 2j * np.pi * freq
    a = lam / C
    g = (1 + tR * p) / (1 + tl * p)
    xi = e / np.sqrt(a)
    b = np.sqrt(lam * C)
    A = xi * np.sqrt(p * g)
    B = b * np.sqrt(p / g)
    S = np.sqrt(sap.lam_z * sap.rho_c) * np.sqrt(p)
    L = R + 1 / S
    if conv == "entropic":
        L = L / (1 + tl * p)
    E = np.exp(-2 * A)
    r = (B * L - 1) / (B * L + 1)
    return (1 + r * E) / ((1 - r * E) * B)


freq = np.logspace(4, np.log10(2e8), 25)
out = {}
for conv in ("local", "entropic"):
    for (tR, tl) in ((25e-12, 10e-12), (25e-12, 25e-12)):
        alpha = 1.0 / 3.0
        # tau_l = (1+alpha) ell^2 C / lam  ->  ell
        ell = np.sqrt(tl * lam / ((1 + alpha) * C_aln))
        m = tw.Material(lam, C_aln, tau_r=tR, ell=ell, alpha=alpha)
        st = tw.Stack((tw.Layer(m, e, 0.0, tw.GKBoundary(temperature=conv)),), sap)
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            zt = tw.kernel(freq, 0.0, st)
        zc = closed(freq, lam, C_aln, e, tR, tl, conv)
        out[f"{conv}:tR={tR},tl={tl}"] = float(np.max(np.abs(zt / zc - 1)))
print(json.dumps(out, indent=1))

# phase difference between GK on the resonance (tR = tl = 25 ps) and the Fourier film, both conventions
for conv in ("local", "entropic"):
    zg = closed(freq, lam, C_aln, e, 25e-12, 25e-12, conv)
    zf = closed(freq, lam, C_aln, e, 0.0, 0.0, conv)
    dphi = np.degrees(np.angle(zg / zf))
    print(conv, "max |dphi| GK-resonance vs Fourier (deg):", float(np.max(np.abs(dphi))))
