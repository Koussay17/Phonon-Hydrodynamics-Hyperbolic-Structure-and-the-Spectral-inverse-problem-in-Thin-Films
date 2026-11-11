(* Workstream B: independent cross-check (Wolfram Language) of the low-frequency Laurent
   coefficients, of the 4x4 moment minors, and of the two T-determinants used in Theorem N2.
   Evaluated on 2026-11-11 in a Wolfram Language kernel; output copied verbatim in
   07w_wolfram_crosscheck.out. *)

ClearAll["Global`*"];
g = (1 + tR p)/(1 + tl p); z = xi^2 p g;
cC = Cosh[Sqrt[z]]; sS = Sinh[Sqrt[z]]/Sqrt[z];
X = xi/b g sS; Y = xi b p sS; S = bs s;
F = (cC + X S)/(Y + cC S) /. p -> s^2;
ser = Series[F, {s, 0, 6}];
cs = Table[Simplify[SeriesCoefficient[ser, k]], {k, -1, 6}];
grads = Table[D[cs[[k + 2]], {{xi, b, tR, tl}}], {k, 0, 6}];
minor[rows_] := Factor[Det[grads[[rows + 1]]]];
{Factor[cs[[1]]], Factor[cs[[2]]], Factor[cs[[3]]], Factor[cs[[4]]], Factor[cs[[5]]],
 minor[{0, 1, 2, 3}], minor[{0, 1, 2, 4}], minor[{0, 2, 3, 4}], Factor[minor[{0, 2, 3, 4}] /. b -> bs]}

(* second evaluation *)
rowsT[ks_, lam_] := Module[{c = Cos[# Pi/4] & /@ ks, sn = Sin[# Pi/4] & /@ ks},
   {c, sn, lam^ks c, lam^ks sn}];
{Factor[Simplify[Det[rowsT[{0, 1, 2, 4}, lam]]]], Factor[Simplify[Det[rowsT[{0, 2, 3, 4}, lam]]]],
 Factor[Simplify[Det[rowsT[{0, 1, 2, 3}, lam]]]]}
