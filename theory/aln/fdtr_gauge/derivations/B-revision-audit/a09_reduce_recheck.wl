(* a09 -- two exact real eliminations re-run on a remote Wolfram Language kernel (12 Nov 2026), inputs copied verbatim
   from derivations/B-revision/b2_02_matrices.wl (b_s = 1).  Outputs in a09_reduce_recheck.out.
   (A) T_eps_Rknown, rank < 4: resolves the Root[...] elided in b2_03_reduce_record.md.
   (B) theta_e read-out, tied eps: the locus rank <= n - 2 (not decided in the record), on the degenerate sets that the
       record decided ({Delta, b = 1, R = 0} for R known; {Pi, R = 0} and {Delta, R = 0, t > 0} for R unknown),
       and the exact null spaces at eps = 1, Fourier corner, b = 1, R = 0. *)

(* ---- (A) *)
M = {{0, 2*tR*tl/b, -tl, tR}, {0, 2*tR/b + 2*tl/b, -1, 1}, {0, 2/b, 0, 0}, {-2*R*b, 0, 0, -eps},
     {eps^2*tR*tl^2, 0, 0, 0}, {R^2*b^2*tl + eps^2*tl^2 + 2*eps*tR*tl, 0, 0, 0},
     {R^2*b^2 - b^2*tl + 2*eps*tl + tR, 0, 0, 0}, {1 - b^2, 0, 0, 0}};
ms = DeleteCases[Union[Expand /@ Flatten[Minors[M, 4]]], 0];
dom = b > 0 && tR >= 0 && tl >= 0 && R >= 0;
res = TimeConstrained[Reduce[And @@ Thread[ms == 0] && dom, {b, tR, tl, R, eps}, Reals], 50, "TIME"];
{Length[ms], res // InputForm, FullSimplify[res, dom] // InputForm}

(* ---- (B) *)
Mk = {{0, 2*tR*tl^2/b, -tl^2, -tR*tl}, {0, 4*tR*tl/b + 2*tl^2/b, -2*tl, -tR - tl}, {0, 2*tR/b + 4*tl/b, -1, -1},
      {0, 2/b, 0, 0}, {-2*R*b*tl, 0, 0, 2*R*b*xi}, {-2*R*b, 0, 0, 1 - eps}, {eps^2*tR*tl^3, 0, 0, -eps^2*tR*tl^2*xi},
      {R^2*b^2*tl^2 + eps^2*tR*tl^2 + eps^2*tl^3 + 2*eps*tR*tl^2, 0, 0, -R^2*b^2*tl*xi - eps^2*tl^2*xi - 2*eps*tR*tl*xi},
      {2*R^2*b^2*tl - b^2*tl^2 + eps^2*tl^2 + 2*eps*tR*tl + 2*eps*tl^2 + tR*tl, 0, 0,
       -R^2*b^2*xi + b^2*tl*xi - 2*eps*tl*xi - tR*xi},
      {R^2*b^2 - 2*b^2*tl + 2*eps*tl + tR + tl, 0, 0, b^2*xi - xi}, {1 - b^2, 0, 0, 0}};
Mu = {{0, 2*tR*tl^2/b, -tl^2, -tR*tl, 0}, {0, 4*tR*tl/b + 2*tl^2/b, -2*tl, -tR - tl, 0},
      {0, 2*tR/b + 4*tl/b, -1, -1, 0}, {0, 2/b, 0, 0, 0}, {-2*R*b*tl, 0, 0, 2*R*b*xi, 0}, {-2*R*b, 0, 0, 1 - eps, 0},
      {eps^2*tR*tl^3, 0, 0, -eps^2*tR*tl^2*xi, 0},
      {R^2*b^2*tl^2 + eps^2*tR*tl^2 + eps^2*tl^3 + 2*eps*tR*tl^2, 0, 0, -R^2*b^2*tl*xi - eps^2*tl^2*xi - 2*eps*tR*tl*xi,
       b*eps*tl^3},
      {2*R^2*b^2*tl - b^2*tl^2 + eps^2*tl^2 + 2*eps*tR*tl + 2*eps*tl^2 + tR*tl, 0, 0,
       -R^2*b^2*xi + b^2*tl*xi - 2*eps*tl*xi - tR*xi, 2*b*eps*tl^2 + b*tl^2},
      {R^2*b^2 - 2*b^2*tl + 2*eps*tl + tR + tl, 0, 0, b^2*xi - xi, b*eps*tl + 2*b*tl}, {1 - b^2, 0, 0, 0, b}};
mk = DeleteCases[Union[Expand /@ Flatten[Minors[Mk /. {xi -> 1, tR -> t, tl -> t, b -> 1, R -> 0}, 3]]], 0];
r1 = Reduce[And @@ Thread[mk == 0] && t >= 0, {t, eps}, Reals];
mu = DeleteCases[Union[Expand /@ Flatten[Minors[Mu /. {xi -> 1, tR -> 0, tl -> 0, R -> 0}, 4]]], 0];
r2 = Reduce[And @@ Thread[mu == 0] && b > 0, {b, eps}, Reals];
mu2 = DeleteCases[Union[Expand /@ Flatten[Minors[Mu /. {xi -> 1, tR -> t, tl -> t, R -> 0}, 4]]], 0];
r3 = Reduce[And @@ Thread[mu2 == 0] && b > 0 && t > 0, {b, t, eps}, Reals];
{r1 // InputForm, r2 // InputForm, r3 // InputForm,
 NullSpace[Mk /. {xi -> 1, tR -> 0, tl -> 0, b -> 1, R -> 0, eps -> 1}],
 NullSpace[Mu /. {xi -> 1, tR -> 0, tl -> 0, b -> 1, R -> 0, eps -> 1}]}
