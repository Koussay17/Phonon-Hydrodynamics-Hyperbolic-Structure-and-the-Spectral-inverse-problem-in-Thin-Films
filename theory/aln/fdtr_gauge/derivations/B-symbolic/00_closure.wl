(* Task 0, independent check in the Wolfram Language.
   Quadrupole of the film closed by the substrate impedance Zs. *)
m = {{Cosh[a], Sinh[a]/bb}, {bb Sinh[a], Cosh[a]}};
front[zs_] := (m[[1, 1]] zs + m[[1, 2]])/(m[[2, 1]] zs + m[[2, 2]]);
zClosed = (bb + s Tanh[a])/(bb (s + bb Tanh[a]));
zBrief = (bb Tanh[a] + s)/(bb (bb + s Tanh[a]));
r1 = FullSimplify[front[1/s] - zClosed];
r2 = FullSimplify[zBrief zClosed bb^2 - 1];
r3 = FullSimplify[front[rr + 1/s] - (bb (rr + 1/s) + Tanh[a])/(bb (1 + bb (rr + 1/s) Tanh[a]))];
thin = {Limit[zClosed, a -> 0], Limit[zBrief, a -> 0]};
(* the film functions of note 7, p = I w, with a numerical spot check against the matrix product *)
vals = {xi1 -> 23/10, b -> 579/100, tR -> 9639/1000, tl -> 6889/100000, bs -> 166/1000, w -> 2663/10000};
aa = xi1 Sqrt[p (1 + tR p)/(1 + tl p)]; bbb = b Sqrt[p (1 + tl p)/(1 + tR p)];
num = N[((Cosh[aa] (1/(bs Sqrt[p])) + Sinh[aa]/bbb)/(bbb Sinh[aa] (1/(bs Sqrt[p])) + Cosh[aa])) /. p -> I w /. vals, 50];
cls = N[(zClosed /. {a -> aa, bb -> bbb, s -> bs Sqrt[p]}) /. p -> I w /. vals, 50];
brf = N[(zBrief /. {a -> aa, bb -> bbb, s -> bs Sqrt[p]}) /. p -> I w /. vals, 50];
{"quad-closed" -> r1, "brief*closed*B^2-1" -> r2, "with R" -> r3, "thin film {closed, brief}" -> thin,
 "spot |quad-closed|/|Z|" -> Abs[num - cls]/Abs[cls], "spot |brief-closed|/|Z|" -> N[Abs[brf - cls]/Abs[cls], 10]}
