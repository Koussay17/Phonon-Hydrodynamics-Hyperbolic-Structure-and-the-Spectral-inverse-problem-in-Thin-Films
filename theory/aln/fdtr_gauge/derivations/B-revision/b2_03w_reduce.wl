(* B-2 / B-4, step 3 (independent system): exact real degenerate sets of the kernel systems of Lemma K.
   Input: b2_02_matrices.wl (b_s = 1). Scale invariance (xi -> mu xi, tau -> mu^2 tau, R -> mu R) allows xi = 1.
   For each model: all maximal minors of Mmat; Reduce over the reals on {all minors = 0} with the domain constraints
   b > 0, tR >= 0, tl >= 0, R >= 0 (eps, epsz, tn real). The output is the semi-algebraic degenerate set
   {theta : rank Mmat(theta) < n}, to be compared with the hand classification of proofs/B-revision.md. *)
SetDirectory[DirectoryName[$InputFileName]];
Get["b2_02_matrices.wl"];
out = OpenWrite["b2_03w_reduce.out"];
wr[s__] := (WriteString[out, s, "\n"]; WriteString[$Output, s, "\n"]);
dom = b > 0 && tR >= 0 && tl >= 0 && R >= 0;
models = Keys[mats];
Do[
  M = mats[name] /. xi -> 1;
  n = Length[First[M]];
  ms = DeleteCases[Union[Expand /@ Flatten[Minors[M, n]]], 0];
  ms = DeleteDuplicates[ms, Expand[#1 + #2] === 0 &];
  vars = Select[{b, tR, tl, R, eps, epsz, tn, to}, ! FreeQ[M, #] &];
  res = TimeConstrained[
    MemoryConstrained[Reduce[And @@ Thread[ms == 0] && dom, vars, Reals], 1500000000, "MEMORY"], 900, "TIME"];
  wr["MODEL ", name, "  n = ", ToString[n], "  minors: ", ToString[Length[ms]], "  vars: ", ToString[vars]];
  wr["  degenerate set: ", ToString[res, InputForm]];
  wr[""],
  {name, models}];
Close[out];
