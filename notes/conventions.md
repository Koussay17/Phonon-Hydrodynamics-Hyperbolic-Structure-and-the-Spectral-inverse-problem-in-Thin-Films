# Conventions retenues

Ce fichier fixe les conventions du projet. Il fait autorité sur toute expression employée dans `src/`.
Toute source extérieure doit être convertie vers ces conventions avant usage.

Dernière mise à jour : 16 septembre 2026, après dérivation de la fermeture grise.

---

## 1. Grandeur de flux

**Convention retenue : densité de flux, par unité de surface, en W·m⁻².**

| Source | Grandeur | Unité | Section apparente |
|---|---|---|---|
| Maillet et al. (2000) | flux total | W | oui, facteur `S` |
| Krapez (2019, en ligne 2018) | densité de flux | W·m⁻² | non |
| Krapez & Rigollet (2017) | densité de flux | W·m⁻² | non |
| Présent travail | densité de flux | W·m⁻² | non |

Motifs : le problème est unidimensionnel, la section n'y joue aucun rôle physique ; et les
comparaisons avec Krapez se font sans conversion.

**Règle d'usage.** Toute expression tirée de Maillet portant un facteur `S` s'emploie avec `S = 1`.

## 2. Vérification effectuée sur le mur homogène

Maillet, §1.3.4 page 10, avec `k = √(p/a)`, `e` l'épaisseur, `S` la section :

```
A = D = cosh(ke)        B = sinh(ke)/(λ k S)        C = λ k S sinh(ke)
```

Expression obtenue par la voie transformée, avec `ξ₁ = e/√a` :

```
A = D = cosh(√p ξ₁)     B = sinh(√p ξ₁)/(b √p)      C = b √p sinh(√p ξ₁)
```

Correspondance, par les identités `λ/√a = b` et `ξ₁ = e/√a` :

```
k e = √(p/a)·e = √p·ξ₁              λ k = (λ/√a)·√p = b √p
```

**Conclusion.** Identiques au seul facteur `S` près. Aucun autre écart : orientation de l'axe, sens du
vecteur d'état et signe du flux coïncident.

## 3. Correspondance avec les notations photothermiques

| Note 2017 | Projet | Grandeur |
|---|---|---|
| `L₂` | `e` | épaisseur du revêtement |
| `α₂` | `a` | diffusivité |
| `κ₂` | `λ` | conductivité |
| `ρ₂C₂` | `ρc` | capacité volumique |
| `P₂` | `b` | effusivité |
| `Q₂` | `ξ₁` | racine du temps de transit, `Q₂ = L₂/√α₂` |
| `b₃₂` | — | rapport d'effusivité substrat sur revêtement |

`Q₂` coïncide exactement avec la coordonnée de Liouville évaluée en face arrière, puisque
`ξ(z) = ∫dz/√a` se réduit à `e/√a` en milieu homogène. **Le couple `{Q₂, P₂}` est le couple
`{ξ₁, b}` du projet**, et le quadripôle de l'équation (1) de la note est celui de la section 2.

## 4. Point de vigilance sur le déterminant

```
A D − B C = cosh² − [sinh/(λkS)]·[λkS sinh] = cosh² − sinh² = 1
```

Le facteur `S` se simplifie. **Le test du déterminant ne détecte donc pas un mélange de conventions.**
Il reste nécessaire mais non suffisant ; la validation exige la comparaison terme à terme.

Seconde limite : sur un empilement épais à haute fréquence, les entrées atteignent 10⁴⁰ et le
déterminant perd toute précision par annulation.

## 5. Impédance de fermeture

**Substrat semi-infini, par unité de surface :**

```
Z = 1 / (b √p)
```

Vérification croisée : Maillet §1.4.2 remarque 2 page 15, et Krapez équation (29).

Réponse en surface d'un système revêtement sur substrat :

```
θ₀ = ℘ (A Z + B) / (C Z + D + h (A Z + B))
```

Cas adiabatique : `h = 0`.

## 6. Coordonnée et grandeurs dérivées

```
ξ(z) = ∫₀ᶻ du/√(a(u))        [ξ] = s^(1/2)
a = λ/ρc                      [a] = m²·s⁻¹
b = √(λ ρc)                   [b] = W·s^(1/2)·m⁻²·K⁻¹
λ = b √a                      ρc = b/√a
```

Forme normale et potentiel :

```
d²ψ/dξ² − (V(ξ) + p) ψ = 0        V = s″/s        s = b^(±1/2)
```

Exposant positif pour la formulation en température, négatif pour la formulation en flux.

## 7. Lois constitutives non-Fourier

Hypothèses : transfert 1D, coefficients constants, perturbation initiale nulle,
pas de source intérieure, conditions de bord du modèle de face avant.

```
Cattaneo : τ_R ∂q/∂t + q = −λ ∂T/∂z
GK historique : τ_R ∂q/∂t + q = −λ ∂T/∂z + 3ℓ² ∂²q/∂z²
λ_eff(p) = λ (1 + τ_ℓ p)/(1 + τ_R p)       τ_ℓ = 3ℓ²/a
σd = ξ₁ √[p(1+τ_R p)/(1+τ_ℓ p)]
λ_eff σ = b √[p(1+τ_ℓ p)/(1+τ_R p)]
```

Remplacer `p` par `p(1+τ_R p)` dans tout le quadripôle de Fourier donne le mauvais
coefficient de flux. La phase du demi-espace de Cattaneo est
`φ = −π/4 + arctan(ωτ_R)/2`, et non une phase tendant vers −π/2.


### 7.1 Convention de fermeture après dérivation

La note 14 établit la fermeture grise conservatrice : le terme 3D est
ell² [laplacien(q) + (1/3) grad(div q)], au lieu du coefficient historique 2.
Poser L²=(1+alpha)ell² et tau_l=L²/a :
historique alpha=2, conservatrice alpha=1/3.
La structure du solveur en tau_l ne change pas. Dans la limite collective,
tau_l=9 tau_N/5 pour la convention historique et 4 tau_N/5 pour la
conservatrice. Ni 1,8 ni 0,8 n'est un rapport hydrodynamique fortement séparé.
Le temps résistif volumique ne doit pas compter à nouveau les frontières déjà
décrites par leurs conditions de bord.
L'API d'entropie ell_sq conserve L²=3 ell_sq : pour une longueur CE
physique, fournir ell_sq=4 ell_CE²/9.

## 8. Invariance et portée

Pour les paramètres homogènes libres, `d→μd, λ→μλ, C→C/μ, ℓ→μℓ`
laisse `b, ξ₁, τ_R, τ_ℓ` et la réponse invariants. Elle ne permet donc pas de
séparer simultanément conductivité, capacité et épaisseur.
Une valeur indépendante de `v` avec la relation imposée `a=v²τ_R/3` peut
restreindre cette invariance. Les profils gradués non-Fourier ne sont pas implémentés.

Une couche graduée Fourier exige les deux propriétés arrière. Si la diffusivité
varie, fournir `Sample.graded_xi1` ou `Layer.graded_xi1` :
les valeurs aux faces ne déterminent pas `∫dz/√a(z)`.
Sans cet argument, l'inférence `d/√a` suppose une diffusivité constante dans la couche,
et est autorisée seulement si les diffusivités aux faces coïncident.
La cohérence avec un profil physique intérieur reste à établir par l'utilisateur.

## 9. Résonance de Fourier

`τ_R=τ_ℓ` redonne le quadripôle de Fourier sous les hypothèses de la section 7.
C'est la **résonance de Fourier** déjà décrite par
[Kovács (2018)](https://arxiv.org/abs/1804.05225).
La reproduction numérique n'établit pas une nouveauté.
[Hennessy et Myers (2021)](https://doi.org/10.1007/978-3-030-64272-3_2)
traite directement GK en thermoréflectance ; le manuscrit accepté est comparé dans la note 15.

Dans la convention historique, avec les coefficients isotropes `a=v²τ_R/3` et `ℓ²=v²τ_Nτ_R/5`
([Lebon et Dauby, 1990](https://doi.org/10.1103/PhysRevA.42.4710)),
`τ_ℓ=9τ_N/5`. Le rapport formel `τ_R/τ_N=1,8` ne satisfait pas `τ_N≪τ_R`.
Il ne prouve pas l'ouverture hydrodynamique lors d'un refroidissement.
Les sources intérieures et conditions initiales générales ne sont pas couvertes
par une affirmation d'indiscernabilité universelle.

## 10. Admissibilité

L'entropie physique `s=s_eq(u)−τ_R q²/(2λT₀²)` est concave.
Son opposé est une entropie mathématique convexe. À l'ordre **quadratique**
en perturbations autour de T₀ :

```
Cattaneo : σ_s = q²/(λT₀²),                 J_s = q/T
GK :       σ_s = [q²+3ℓ²(∂q/∂z)²]/(λT₀²), J_s = q/T+3ℓ²q(∂q/∂z)/(λT₀²)
```

Garder `T=T₀+εθ` et `q=εj` avant le développement ; remplacer T par une
constante trop tôt efface les gradients et peut masquer une erreur.
La production est non négative pour `λ>0, ℓ²≥0`; la concavité stricte en flux
exige `τ_R>0`.

| Loi | Propagation dans le modèle continu |
|---|---|
| Fourier | infinie |
| Cattaneo, τ_R>0 | finie, √(a/τ_R) |
| GK, τ_R>0, ℓ>0 | infinie ; ℓ=0 redonne Cattaneo |

Avec `exp(iωt+ikz)`, choisir `k=iσ` pour une décroissance vers z>0.
Le système principal de Cattaneo en (u,q) a la matrice
`[[0,1],[a/τ_R,0]]`. Les vitesses de phase seules ne remplacent pas
l'analyse du support causal des solutions.

## 11. Données matériau et limites

### 11.1 Conductivité et contraste

[Cheng et al. (2020)](https://arxiv.org/abs/1911.01595) mesure environ
321 W·m⁻¹·K⁻¹ sur des films MOCVD de **18 et 22,5 µm**.
Cette valeur ne doit pas être attribuée sans mesure à un film de 500 nm.
Le défaut du code, 60 W·m⁻¹·K⁻¹, reste un scénario illustratif.

Pour le substrat supposé λ_s=35, C_s=3,03×10⁶, le coefficient
Fourier/perfect-contact `Γ=(b_f−b_s)/(b_f+b_s)` vaut 0,4596 pour λ_f=321,
environ 0,077 pour λ_f=60, et zéro à λ_f≈44 lorsque C_f=2,41×10⁶.
Le risque pratique dépend des propriétés du film réel et des contacts.

### 11.2 Temps modaux

[Ma, Li et Luo (2014)](https://doi.org/10.1103/PhysRevB.90.035203), figure 7 :
lectures approximatives (facteur d'environ deux), sur une portion acoustique
de l'axe de **pulsation** ω étiqueté THz, de 2 à 10.

| Branche | ω (axe source) | τ_U (ps) | τ_N (ps) | τ_U/τ_N |
|---|---|---|---|---|
| LA | 2 | 7×10⁵ | 10⁴ | 70 |
| LA | 10 | 10⁴ | 4×10² | 25 |
| TA | 2 | 4×10⁵ | 2×10³ | 200 |
| TA | 10 | 4×10³ | 4×10² | 10 |

Ce ne sont ni des données numérisées précises, ni tout le spectre, ni des temps
effectifs GK. L'accord de pentes n'auto-valide pas les lectures.
τ_U≈τ_R est une approximation ici ; les autres processus résistifs sont omis.
À basse pulsation, Ma et al. donne des taux en ω³ pour U, ω pour N-TA et ω² pour N-LA.

### 11.3 Carte des hiérarchies

`x=τ_R/τ_N`, `y=d/(vτ_N)`. La fenêtre forte exige `1≪y≪x`.
`y<1` seul ne suffit pas : les frontières dominent si `y<min(1,x)`.
Le script 05 conserve les couples (x,y) du même mode.
À 500 µm, le point TA d'abscisse 10 a y≈208 : tous les points ne sont pas
hydrodynamiques. À 500 nm, y≈0,0083–0,208 pour ces lectures seulement.
Aucune trajectoire en température n'est établie.

### 11.4 Longueur grise et longueur normale

vτ_N≈2,4–60 µm pour ces points et v=6000 m/s.
`Λ_g=3λ/(Cv)` est une estimation grise liée à la conductivité, pas le libre
parcours de toutes les collisions ni une reconstruction du spectre.

| Hypothèses à d=500 nm, C=2,41×10⁶, v=6000 | Λ_g | Kn_g |
|---|---|---|
| λ=60 | 12,45 nm | 0,0249 |
| λ=321 | 66,60 nm | 0,1332 |

`mean_free_path` conserve son nom d'API, avec cette définition grise ;
`transport_regime` donne une étiquette heuristique (seuils 0,1 et 1).
Le diagnostic ne valide pas automatiquement la diffusion ou une classification
balistique d'un film réel. Comparer une conductivité apparente au massif est
pertinent pour étudier sa suppression.
[Hoque et al. (2024)](https://arxiv.org/abs/2409.14328) étudie une transition
propre à ses films ; aucun seuil universel de 67 nm n'en découle.

### 11.5 Erreur de fermeture

L'écart Callaway hors plan d'environ −12 % chez Ma et al. concerne la conductivité
stationnaire. Il n'est pas une borne sur l'erreur d'un temps dynamique ajusté.
L'anisotropie d'environ 5 % est celle de la référence ab initio.
Le modèle de Cheng sans processus normaux ne fournit pas τ_N.

## 12. Conception instrumentale et numérique

### 12.1 Fréquence caractéristique

`f_c=1/(2πξ₁²)`. Pour d=500 nm et C=2,41×10⁶ :
15,85 MHz avec λ=60, **84,80 MHz** avec λ=321.
Sous f_c, le substrat peut dominer sans annuler exactement la sensibilité au film.

### 12.2 Interface aveugle

À contact parfait et sous Fourier, les effusivités égales rendent l'interface
invisible. Une résistance de contact ou des temps non-Fourier différents
peuvent lever cette égalité des impédances ; Γ seul n'est pas un diagnostic complet.

### 12.3 Centre de transition, pas seuil de mesurabilité

`∂φ/∂lnτ = ωτ/[2(1+(ωτ)²)]` atteint 1/4 rad à ωτ=1.
À 200 MHz, le temps centré en haut de bande est 0,7958 ns, **pas τ_min mesurable**.
Exemple idéal, τ=10 ps : écart de phase 0,35998° à 200 MHz ;
bruit 0,01°, τ seul inconnu → incertitude relative locale 2,78 %.
Cela ne démontre pas la faisabilité sur un échantillon réel.
La relation φ=arg(E)/2 signifie que les deux angles ne sont pas des
validations indépendantes.

### 12.4 Corrélations et covariance

Libérer un temps peut augmenter les incertitudes des autres paramètres.
La note IV donne des facteurs 3,13 et 1,61 sur λ et C.
Les séparations de 16 % et 25 % autrefois annoncées sont retirées.
Le script 04 utilise **10 kHz–200 MHz, 80 points, 1 % d'amplitude, 0,1° de phase**,
avec quatre paramètres libres.

`fisher_analysis` et sa variante temporelle lèvent `NonIdentifiableError`
si le rang ne permet pas une covariance séparée. La SVD des colonnes normalisées
distingue une faible sensibilité d'une dépendance. Les ajustements peuvent
converger avec une covariance infinie : `success` décrit l'optimiseur, pas
l'identifiabilité. Le conditionnement utilise les paramètres logarithmiques.

### 12.5 Débordement

La croissance de cosh/sinh dépend de **|Re(σd)|**, pas de |σd|.
Les réponses homogènes et les empilements Fourier homogènes utilisent une
récursion d'impédance avec tanh saturée, stable pour les couches très épaisses.
Les matrices explicites restent limitées ; les couches graduées très épaisses
ne bénéficient pas encore de cette stabilisation complète.
Les matrices ont leur limite analytique à p=0 ; l'impédance du demi-espace
à fréquence nulle est non bornée.

## 13. Inversion de Laplace

Gaver–Stehfest est implémenté, n=12 par défaut. Les erreurs rapportées sur deux
décades concernaient un benchmark exponentiel, pas une garantie pour tout signal.
Vérifier le signal **et ses sensibilités**, particulièrement près d'un front d'onde.
De Hoog est une piste de remplacement, **pas un repli déjà implémenté**.

## 14. Validation

Les essais synthétiques utilisant le même modèle direct vérifient l'estimateur,
sans valider une expérience. L'accord avec Monte Carlo est conditionnel au
modèle et au bruit. Composition, déterminant et limites ne suffisent pas.

Ancrages : Maillet et al., mur Fourier ; réponse impulsionnelle `Q/(b√(πt))` ;
phase du demi-espace Cattaneo de Camacho de la Rosa et al. (2025).
L'intégration indépendante des bilans différentiels vérifie aussi les quadripôles
homogènes GK et gradués Fourier sans réutiliser leurs formules fermées.

## 15. Reproductibilité et suite

Installer `requirements-dev.txt`, lancer pytest puis `scripts/rebuild.py`.
`seed` contrôle l'initialisation et emcee. `autocorr_reliable` indique si la
chaîne suffit à l'estimation d'autocorrélation ; il ne certifie pas toute convergence.
Les priors log-uniformes bornés sont des hypothèses explicites.

La fermeture grise est dérivée dans la note 14, les comparaisons sont dans la
note 15, l'étude expérimentale synthétique dans la note 16.
La note 17 ajoute le spectre RTA à 300 K et une trajectoire massive publiée.
Restent : fermeture dynamique de l'AlN et trajectoire N/U en température,
géométrie/calibration du banc et données réelles. Le texte intégral IJHMT de
Krapez 2016 et le supplément Camacho n'ont pas été obtenus.

## 16. Base FDTR axisymétrique

src/fdtr.py implémente Fourier uniquement, avec absorption surfacique,
faisceaux gaussiens coaxiaux, couches homogènes éventuellement anisotropes et
substrat semi-infini. Rayons à 1/e² ; puissance absorbée totale ; réponse en K/W.
La phase est l'argument de la moyenne complexe, pas une moyenne de phases.
Cette base ne contient pas les conditions de bord hydrodynamiques 3D de Beardo.
scripts/analyse_experiment.py produit uniquement des résultats synthétiques
définis dans theory/experimental_results.json et la figure 06.

## 17. Traitement spectral AlN

Les données Rao comportent 12 branches à 300 K. Les taux anharmoniques
ne séparent pas normal et umklapp ; plus/minus désignent absorption/émission.
Les temps statique et mémoire ont des pondérations différentes et ne sont
pas des temps GK identifiés. Les opérateurs conservateurs demandent la zone
complète ou une reconstruction explicite des symétries.

La loi de suppression de src/spectral.py concerne la conduction parallèle
entre surfaces réfléchissantes, sous RTA stationnaire. Pas de double comptage
des frontières, ni de transposition directe à FDTR transversal.
La courbe kappa(T) publiée et le spectre Rao sont deux calculs indépendants.
C(T) est harmonique, à fréquences fixes ; les durées ne sont pas extrapolées.
