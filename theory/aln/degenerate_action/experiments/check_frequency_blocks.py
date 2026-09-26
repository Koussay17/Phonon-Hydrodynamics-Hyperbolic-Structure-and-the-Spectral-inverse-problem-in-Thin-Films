"""Regression checks for the independently discovered overlap defect."""
from pathlib import Path
import json,numpy as np
from frequency_blocks import frequency_blocks
def groups(x,t=1e-10):
 return [g.tolist() for g in frequency_blocks(x,t)]
def rejects(x,t):
 try: frequency_blocks(x,t)
 except ValueError: return True
 return False
checks={
 "exact_doublet":groups([1.,1.,2.])==[[0,1],[2]],
 "unsorted_partition":groups([2.,1.,1.])==[[1,2],[0]],
 "chain_rejected":rejects([1.,1.75,2.5],1.),
 "long_chain_rejected":rejects([1.,1.75,2.5,3.25],1.),
 "threshold_boundary_separate":groups([1.,2.],1.)==[[0],[1]],
 "nonfinite_rejected":rejects([1.,float("nan")],1.),
 "invalid_tolerance_rejected":rejects([1.],0.),
 "empty_partition":groups([])==[],
}
assert all(checks.values()),checks
print(json.dumps({"checks":checks},indent=2))
