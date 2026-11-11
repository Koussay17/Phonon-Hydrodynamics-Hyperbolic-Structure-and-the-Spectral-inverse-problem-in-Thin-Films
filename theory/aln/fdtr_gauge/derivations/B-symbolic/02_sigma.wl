(* Task 2, independent check in the Wolfram Language.
   Z from the quadrupole product; identities of the Jacobian columns; determinant factorisation;
   ranks on the strata (singular values without equilibration, 50 digits). *)
z[p_, x1_, b_, tR_, tl_, bs_] := Module[{aa, bb, zs, m},
  aa = x1 Sqrt[p (1 + tR p)/(1 + tl p)]; bb = b Sqrt[p (1 + tl p)/(1 + tR p)];
  zs = 1/(bs Sqrt[p]); m = {{Cosh[aa], Sinh[aa]/bb}, {bb Sinh[aa], Cosh[aa]}};
  (m[[1, 1]] zs + m[[1, 2]])/(m[[2, 1]] zs + m[[2, 2]])];
zz = z[p, x1, b, tR, tl, bs];
{d1, d2, d3, d4} = D[zz, {{x1, b, tR, tl}}];
gR = p/(2 (1 + tR p)); gl = p/(2 (1 + tl p));
asm = p > 0 && x1 > 0 && b > 0 && bs > 0 && tR >= 0 && tl >= 0;
i1 = Simplify[d3 - gR (x1 d1 - b d2), asm];
i2 = Simplify[d4 + gl (x1 d1 - b d2), asm];
i3 = Simplify[(d3 + d4) /. tl -> tR, asm];
i4 = Simplify[d1 /. {tl -> tR, b -> bs}, asm];
(* the same identities at p = I w, 50 digits, random rational points *)
SeedRandom[5];
num = Table[Module[{r = {x1 -> RandomInteger[{1, 300}]/100, b -> RandomInteger[{1, 300}]/100,
      tR -> RandomInteger[{0, 300}]/100, tl -> RandomInteger[{0, 300}]/100, bs -> RandomInteger[{1, 300}]/100,
      p -> I RandomInteger[{1, 3000}]/100}},
    N[{Abs[d3 - gR (x1 d1 - b d2)]/Abs[d3], Abs[d4 + gl (x1 d1 - b d2)]/Abs[d4]} /. r, 50]], {10}];
(* determinant factorisation with free U_m, V_m, p_m *)
row[m_] := {u[m]/x1, v[m]/b, (p[m]/(2 (1 + tR p[m]))) (u[m] - v[m]), -(p[m]/(2 (1 + tl p[m]))) (u[m] - v[m])};
lhs = Det[Table[row[m], {m, 4}]];
kap[q_] := q^2/(2 (1 + tR q) (1 + tl q));
rhs = -(tR - tl)/(x1 b) Product[u[m] - v[m], {m, 4}] Det[Table[{1, v[m]/(u[m] - v[m]), p[m]/(2 (1 + tR p[m])), kap[p[m]]}, {m, 4}]];
fac = Together[lhs - rhs];
(* ranks on strata: singular values of the real Jacobian, 6 frequencies, no equilibration *)
jr[pt_, oms_] := Module[{rows = Table[N[{d1, d2, d3, d4} /. pt /. p -> I om, 50], {om, oms}]},
   Flatten[{Re[#], Im[#]} & /@ rows, 1]];
oms = {1/10, 1/2, 1, 3, 10, 40};
svn[mat_] := Module[{s = SingularValueList[mat, Tolerance -> 0]}, N[s/First[s], 5]];
strata = {
  "generic" -> svn[jr[{x1 -> 1, b -> 3, tR -> 1/2, tl -> 1/5, bs -> 1}, oms]],
  "tau_R = tau_l, b != bs" -> svn[jr[{x1 -> 1, b -> 3, tR -> 1/2, tl -> 1/2, bs -> 1}, oms]],
  "tau_R = tau_l, b = bs" -> svn[jr[{x1 -> 1, b -> 3, tR -> 1/2, tl -> 1/2, bs -> 3}, oms]],
  "Fourier tau = 0, b != bs" -> svn[jr[{x1 -> 1, b -> 3, tR -> 0, tl -> 0, bs -> 1}, oms]],
  "Fourier tau = 0, b = bs" -> svn[jr[{x1 -> 1, b -> 3, tR -> 0, tl -> 0, bs -> 3}, oms]],
  "b = bs, tau_R != tau_l" -> svn[jr[{x1 -> 1, b -> 3, tR -> 1/2, tl -> 1/5, bs -> 3}, oms]],
  "Cattaneo tau_l = 0" -> svn[jr[{x1 -> 1, b -> 3, tR -> 1/2, tl -> 0, bs -> 1}, oms]],
  "tau_R = 0" -> svn[jr[{x1 -> 1, b -> 3, tR -> 0, tl -> 1/2, bs -> 1}, oms]]};
nulls = {"null diag" -> Chop[N[(jr[{x1 -> 1, b -> 3, tR -> 1/2, tl -> 1/2, bs -> 1}, oms]) . {0, 0, 1, 1}, 30], 10^-40] === Table[0, {12}],
  "null diag b=bs" -> Chop[N[(jr[{x1 -> 1, b -> 3, tR -> 1/2, tl -> 1/2, bs -> 3}, oms]) . {1, 0, 0, 0}, 30], 10^-40] === Table[0, {12}]};
{"i1 (tauR column)" -> i1, "i2 (taul column)" -> i2, "i3 (diag sum)" -> i3, "i4 (dZ/dxi1 at b=bs, diag)" -> i4,
 "max numeric residuals" -> Max[num], "det factorisation" -> fac, "strata" -> strata, "nulls" -> nulls}
