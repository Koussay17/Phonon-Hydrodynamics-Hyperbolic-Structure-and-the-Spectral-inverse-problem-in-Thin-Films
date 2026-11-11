(* Task 3, independent check in the Wolfram Language. *)
z[p_, x1_, b_, tR_, tl_, bs_] := Module[{aa, bb, zs, m},
  aa = x1 Sqrt[p (1 + tR p)/(1 + tl p)]; bb = b Sqrt[p (1 + tl p)/(1 + tR p)];
  zs = 1/(bs Sqrt[p]); m = {{Cosh[aa], Sinh[aa]/bb}, {bb Sinh[aa], Cosh[aa]}};
  (m[[1, 1]] zs + m[[1, 2]])/(m[[2, 1]] zs + m[[2, 2]])];
zq = z[p, e/Sqrt[a], lam/Sqrt[a], tR, k ell^2/a, bs];
euler = e D[zq, e] + 2 a D[zq, a] + lam D[zq, lam] + ell D[zq, ell];
asm = p > 0 && e > 0 && a > 0 && lam > 0 && ell > 0 && tR > 0 && k > 0 && bs > 0;
r1 = Simplify[euler, asm];
noell = Simplify[euler - ell D[zq, ell], asm] === 0;
spot = N[Abs[euler]/Abs[e D[zq, e]] /. {p -> I 37/10, e -> 7/10, a -> 13/10, lam -> 21/10, ell -> 3/10, tR -> 9/10, k -> 3, bs -> 11/10}, 40];
spotNoEll = N[Abs[euler - ell D[zq, ell]]/Abs[e D[zq, e]] /. {p -> I 37/10, e -> 7/10, a -> 13/10, lam -> 21/10, ell -> 3/10, tR -> 9/10, k -> 3, bs -> 11/10}, 20];
finite = Simplify[(zq /. {e -> mu e, a -> mu^2 a, lam -> mu lam, ell -> mu ell}) - zq, asm && mu > 0];
phi = {e/Sqrt[a], lam/Sqrt[a], tR, k ell^2/a};
jphi = D[phi, {{e, a, lam, ell, tR}}];
(* dimension exponents (m, kg, s, K), written by hand *)
dims = <|"e" -> {1, 0, 0, 0}, "a" -> {2, 0, -1, 0}, "lambda" -> {1, 1, -3, -1}, "rhoc" -> {-1, 1, -2, -1},
   "ell" -> {1, 0, 0, 0}, "tauR" -> {0, 0, 1, 0}, "bs" -> {0, 1, -5/2, -1}, "omega" -> {0, 0, -1, 0}, "Z" -> {0, -1, 3, 1}|>;
kerData = NullSpace[{dims["omega"], dims["bs"], dims["Z"]}];
action = (Values[KeyTake[dims, {"e", "a", "lambda", "rhoc", "ell", "tauR"}]] . #) & /@ kerData;
{"Euler D Z" -> r1, "Euler without ell term is zero?" -> noell, "spot |DZ|/|e dZ/de|" -> spot,
 "spot without ell" -> spotNoEll, "finite mu" -> finite,
 "rank dphi" -> MatrixRank[jphi], "kernel dphi" -> Simplify[NullSpace[jphi]],
 "unit group preserving (omega, bs, Z)" -> kerData, "its action on ln(e,a,lambda,rhoc,ell,tauR)" -> action}
