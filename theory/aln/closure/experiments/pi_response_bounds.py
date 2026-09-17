"""Independent PI check of real-Laplace scalar response bounds, synthetic SPD only."""
import json
from pathlib import Path
import numpy as np
rng=np.random.default_rng(16092026)
worst=0.; minimum_margin=1.
for case in range(100):
    n=8
    q,_=np.linalg.qr(rng.normal(size=(n,n)))
    rates=np.exp(rng.uniform(-3,3,n))
    c=(q*rates)@q.T
    b=rng.normal(size=n)
    x=np.linalg.solve(c,b)
    k0=b@x
    ts=k0/(b@b)
    tm=(x@x)/k0
    for s in np.geomspace(1e-5,1e3,40):
        exact=b@np.linalg.solve(c+s*np.eye(n),b)
        lower=k0/(1+s*tm)
        upper=k0/(1+s*ts)
        violation=max(lower-exact,exact-upper,0)/k0
        worst=max(worst,violation)
        minimum_margin=min(minimum_margin,(exact-lower)/k0,(upper-exact)/k0)
assert worst<1e-12
out={'seed':16092026,'synthetic_SPD_operators':100,'positive_Laplace_points_each':40,
     'maximum_relative_bound_violation':worst,'minimum_relative_margin':minimum_margin,
     'scope':'positive real Laplace parameter; finite dissipative SPD space; no AlN inference'}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
