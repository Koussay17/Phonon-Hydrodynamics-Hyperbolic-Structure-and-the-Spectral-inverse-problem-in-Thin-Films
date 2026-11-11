(* Task 4, independent check in the Wolfram Language.
   Small p: Series in q = Sqrt[p]; large p: Series of 1/B in u = p^(-1/2); Cattaneo and tau_R = 0
   strata.  Exact values of the coefficients at the rational test point are returned for comparison
   with the SymPy output (04_expansions.py). *)
pt = {x1 -> 7/10, b -> 23/10, bs -> 11/10, tR -> 3/10, tl -> 1/8, r -> 1/5};
zfun[aa_, bb_, zs_, t_] := (bb zs + t)/(bb (1 + bb zs t));
rho = Sqrt[(1 + tR q^2)/(1 + tl q^2)];
zq0 = zfun[x1 q rho, b q/rho, 1/(bs q), Tanh[x1 q rho]];
zqR = zfun[x1 q rho, b q/rho, r + 1/(bs q), Tanh[x1 q rho]];
s0 = Series[zq0, {q, 0, 5}]; sR = Series[zqR, {q, 0, 5}];
c0 = Table[Factor[Simplify[SeriesCoefficient[s0, n]]], {n, -1, 5}];
cR = Table[Factor[Simplify[SeriesCoefficient[sR, n]]], {n, -1, 5}];
(* large p generic: 1/B with p = 1/u^2 *)
invB = (u/b) Sqrt[(1 + tR/u^2)/(1 + tl/u^2)];
sB = Series[invB, {u, 0, 8}, Assumptions -> u > 0 && tR > 0 && tl > 0];
yB = Table[Factor[Simplify[SeriesCoefficient[sB, 2 k + 1], Assumptions -> tR > 0 && tl > 0]], {k, 0, 3}];
(* Cattaneo tau_l = 0: B = b Sqrt[p/(1 + tR p)], S = bs Sqrt[p], t fixed *)
bC = b/Sqrt[tR + u^2]; sC = Series[zfun[0, bC, u/bs, t], {u, 0, 3}];
cC = Table[Factor[Simplify[SeriesCoefficient[sC, k]]], {k, 0, 3}];
(* tau_R = 0 < tau_l *)
bN = b Sqrt[1 + tl/u^2]/u; aN = x1/Sqrt[tl + u^2];
sN = Series[zfun[aN, bN, u/bs, Tanh[aN]], {u, 0, 5}, Assumptions -> u > 0 && tl > 0];
cN = Table[Factor[Simplify[SeriesCoefficient[sN, k]]], {k, 2, 5}];
{"small p, R = 0" -> c0, "values" -> (c0 /. pt),
 "small p, with R" -> cR, "values R" -> (cR /. pt),
 "1/B coefficients of p^(-1/2-k)" -> yB, "values y" -> (yB /. pt),
 "Cattaneo coefficients of p^(-k/2), t fixed" -> cC,
 "tau_R = 0 coefficients of p^(-k/2), k = 2..5" -> cN, "values N" -> N[(cN /. tR -> 0 /. pt), 30]}
