"""Audit helper: closure of an exported event set under the 24 symmetry operations (diagnoses missing images)."""
import sys
import numpy as np
Z = np.load(sys.argv[1])
f = Z["freqs"]; Nq, nb = f.shape; M = Nq * nb
ev_p, ev_a, ev_b, g, dl = Z["ev_p"], Z["ev_a"], Z["ev_b"], Z["ev_g"], Z["ev_dl"]
m = len(g); sig = float(Z["sigma"]); nsig = float(Z["nsig"])
key = (ev_p * M + ev_a) * M + ev_b
order = np.argsort(key); skey = key[order]
missing = np.zeros(m, dtype=bool)
for Rm in Z["rot_maps"]:
    p2 = Rm[ev_p // nb] * nb + ev_p % nb; a2 = Rm[ev_a // nb] * nb + ev_a % nb; b2 = Rm[ev_b // nb] * nb + ev_b % nb
    k2 = (p2 * M + np.minimum(a2, b2)) * M + np.maximum(a2, b2)
    pos = np.minimum(np.searchsorted(skey, k2), m - 1)
    missing |= skey[pos] != k2
print("events", m, "with a missing image", int(missing.sum()))
if missing.any():
    mm = np.nonzero(missing)[0]
    print("|dl|/sigma of those: min %.12f max %.12f" % (np.abs(dl[mm]).min() / sig, np.abs(dl[mm]).max() / sig))
    print("rate of those / median rate: max %.3e" % (g[mm].max() / np.median(g)))
    print("share of total rate: %.3e" % (g[mm].sum() / g.sum()))
    print("fraction within 1e-9 THz of the cut:", float(np.mean(np.abs(np.abs(dl[mm]) - nsig * sig) < 1e-9)))
