"""Disjoint numerical frequency blocks; reject ambiguous chains.

This is a numerical partition rule, not a physical degeneracy certificate.
"""
import numpy as np

def frequency_blocks(frequencies, tolerance=1e-10):
    f=np.asarray(frequencies,dtype=float)
    if f.ndim!=1 or not np.all(np.isfinite(f)):
        raise ValueError("finite one-dimensional frequencies required")
    if not np.isfinite(tolerance) or tolerance<=0:
        raise ValueError("finite positive tolerance required")
    if not len(f):
        return []
    order=np.argsort(f,kind="stable")
    groups=[]
    start=0
    for j in range(1,len(order)+1):
        if j==len(order) or f[order[j]]-f[order[j-1]]>=tolerance:
            group=order[start:j]
            if f[group[-1]]-f[group[0]]>=tolerance:
                raise ValueError("ambiguous chained near-degeneracy; specify physical blocks")
            groups.append(np.sort(group))
            start=j
    return groups
