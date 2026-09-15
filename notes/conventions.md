# Conventions retenues

Ce fichier fixe les conventions du projet. Il fait autorité sur toute expression employée dans `src/`.
Toute source extérieure doit être convertie vers ces conventions avant usage.

Dernière mise à jour : 13 septembre 2026, sixième révision.

---

## 1. Grandeur de flux

**Convention retenue : densité de flux, par unité de surface, en W·m⁻².**

| Source | Grandeur | Unité | Section apparente |
|---|---|---|---|
| Maillet et al. (2000) | flux total | W | oui, facteur `S` |
| Krapez (2018) | densité de flux | W·m⁻² | non |
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

### 7.1 Cattaneo

```
τ_R ∂φ/∂t + φ = −λ ∂T/∂z        d'où        λ_eff = λ / (1 + τ_R p)
```

### 7.2 Guyer–Krumhansl

En une dimension, les deux termes non locaux se combinent en `3 ∂²q/∂z²`, et l'élimination du flux
par le bilan d'énergie donne une **conductivité effective dépendant de la fréquence** :

```
λ_eff(p) = ( λ + 3 ℓ² p ρc ) / ( 1 + τ_R p )
```

### 7.3 Les deux seules combinaisons qui entrent dans la réponse

```
σ·e      = ξ₁ · √[ p (1 + τ_R p) / (1 + τ_ℓ p) ]
λ_eff·σ  = b  · √[ p (1 + τ_ℓ p) / (1 + τ_R p) ]
```

avec

```
τ_ℓ = 3 ℓ² / a          le temps de diffusion sur la longueur non locale
```

Les deux temps **échangent leurs rôles** entre les deux combinaisons. Cattaneo correspond à
`τ_ℓ = 0`, Fourier à `τ_R = τ_ℓ = 0`.

### 7.4 Erreur à ne pas reproduire

**Substituer `p → p(1+τ_R p)` dans l'ensemble du quadripôle de Fourier est faux.**

Cette substitution donne le bon argument des fonctions hyperboliques, mais un coefficient de flux
erroné d'un facteur `(1 + τ_R p)`. La conséquence est visible : la phase d'un milieu semi-infini
part alors de −45° vers −90°, alors que la forme fermée de Camacho de la Rosa et al. (2025) donne
−45° vers 0°, en passant par −22,5° à `ωτ = 1`.

Cette erreur a survécu à cent treize tests de cohérence interne — déterminant unimodulaire, limite
homogène, composition, découpage — parce qu'elle était cohérente avec elle-même. Seule la
confrontation à une forme fermée extérieure l'a révélée.

Elle n'a pas affecté les conclusions d'identifiabilité, le facteur fautif ne dépendant ni de `λ` ni
de `ρc`, donc étant invariant sous le groupe. Elle affectait en revanche toute comparaison à une
mesure réelle.

### 7.5 Lien microscopique

```
ℓ² = v² τ_N τ_R / 5          a = v² τ_R / 3          d'où          τ_ℓ = 9 τ_N / 5
```

## 8. Paramétrage naturel

**Ce que la mesure de face avant détermine est `{b, ξ₁}`, pas `{λ, ρc}`.**

Les vecteurs propres de la matrice de Fisher le confirment : très au-dessus de la transition, la
direction bien déterminée est `(1, 1)` en coordonnées logarithmiques, soit `log λ + log ρc = 2 log b`,
l'effusivité. La direction orthogonale, la diffusivité, porte une valeur propre inférieure de onze
ordres de grandeur.

**Règle.** Estimer en espace logarithmique. Pour les problèmes mal conditionnés, borner l'estimation.

## 9. Identifiabilité — résultats acquis

### 9.1 Sous Fourier

**Non-identifiabilité inconditionnelle** du triplet `{e, a, λ}` en configuration revêtement sur
substrat, excitation et détection en face avant, transfert unidimensionnel. Krapez & Rigollet (2017).

Invariance du groupe à un paramètre

```
e → μe        a → μ²a        λ → μλ        ρc → μ⁻¹ρc
```

sous lequel `ξ₁` et `b` sont invariants, donc le quadripôle aussi. Forme infinitésimale, l'identité
d'Euler :

```
e ∂G/∂e + 2a ∂G/∂a + λ ∂G/∂λ = 0        pour toute fréquence
```

Le cas gradué (Krapez 2023) présente la même invariance, de dimension fonctionnelle au lieu de un :
tout reparamétrage croissant de la profondeur préservant `b(ξ)` laisse la réponse inchangée.

**Levée.** Seule une information imposant une échelle de longueur indépendante la restreint : mesure
directe d'une propriété, épaisseur mesurée indépendamment, diffusion latérale bidimensionnelle,
sources internes.

### 9.2 Sous Cattaneo et sous Guyer–Krumhansl

**L'invariance survit intacte dans les deux cas.**

La raison est structurelle : la réponse ne dépend du film que par `ξ₁`, `b`, `τ_R` et `τ_ℓ`.
**Aucune longueur n'y figure séparément.** La longueur non locale n'entre que par `ℓ²/a`, qui est un
temps, et les temps sont invariants par changement d'échelle en profondeur.

La relation d'Euler doit alors être étendue :

```
e ∂G/∂e + 2a ∂G/∂a + λ ∂G/∂λ + ℓ ∂G/∂ℓ = 0
```

Le terme en `ℓ` est indispensable : résidu de 10⁻⁵ avec, de 0,3 sans. La longueur non locale se
transforme donc bien comme une longueur, `ℓ → μℓ`, ce qui est cohérent avec `v → μv` puisque
`v² = 3a/τ_R`.

Voir `tests/test_cattaneo.py` et `tests/test_guyer_krumhansl.py`.

Les deux temps entrent en revanche comme des paramètres authentiquement indépendants : leur
estimation ne dégrade ni celle de `b` ni celle de `ξ₁`.

### 9.3 Sensibilité et identifiabilité

| | Non-identifiabilité structurelle | Mauvais conditionnement |
|---|---|---|
| Origine | rang déficient | valeurs propres faibles mais non nulles |
| Dépendance au bruit | aucune | déterminante |
| Dépendance à la bande | aucune | déterminante |
| Effet | solution non unique | solution unique, incertitude grande |
| Remède | information de nature différente | élargir la bande, réduire le bruit |

Une dérivée du signal grande n'établit pas qu'un paramètre soit identifiable. La grandeur
déterminante est la composante de la colonne jacobienne **orthogonale** aux autres colonnes, donc un
angle et non une norme :

```
σ(ln θᵢ) = σ_bruit / ‖J⊥,ᵢ‖
```

Un point où toutes les sensibilités croissent simultanément — typiquement une quasi-résonance — est
le plus favorable en norme et le plus défavorable en angle.

## 10. Admissibilité des lois constitutives

### 10.1 Les deux critères sont indépendants

| Loi | Entropie convexe | Seconde loi | Propagation |
|---|---|---|---|
| Fourier | oui | oui, si `λ > 0` | **infinie** |
| Cattaneo | oui, si `τ_R > 0` | oui, si `λ > 0` | **finie**, `√(a/τ_R)` |
| Guyer–Krumhansl | oui, si `τ_R > 0` | oui, si `λ > 0` et `ℓ² ≥ 0` | **infinie** |

Les trois lois sont thermodynamiquement admissibles. **Une seule propage à vitesse finie.** La
hiérarchie Fourier → Cattaneo → Guyer–Krumhansl n'est donc pas un raffinement monotone : le modèle
le plus riche perd la propriété qui motive habituellement l'abandon de Fourier.

### 10.2 Dispersion

```
σ² = p (1 + τ_R p) / [ a (1 + τ_ℓ p) ]          σ = i k
```

Comportement du nombre d'onde à haute fréquence, exposant mesuré :

| Loi | `k ∝ ω^n` | Vitesse de phase |
|---|---|---|
| Fourier | `n = 0,5000` | non bornée |
| Cattaneo | `n = 1,0000` | sature à `√(a/τ_R)` |
| Guyer–Krumhansl | `n = 0,5004` | non bornée |

Le terme non local est **diffusif sur le flux lui-même**. Au-dessus de la fréquence où il domine le
terme en `τ_R ∂q/∂t`, la partie principale du système redevient parabolique.

En modèle de Debye, `a = v²τ_R/3`, donc la vitesse limite de Cattaneo vaut `v/√3`.

### 10.3 Entropie étendue

La thermodynamique irréversible étendue prend le flux comme variable d'état :

```
s(u, q) = s_éq(u) − ( τ_R / 2λT₀² ) · q²
```

concave si `τ_R > 0`. Production, avec le flux d'entropie `q/T` :

```
Cattaneo :          σ_s = q² / (λT₀²)
Guyer–Krumhansl :   σ_s = [ q² + 3ℓ² (∂q/∂z)² ] / (λT₀²)
```

La seconde exige un **flux d'entropie étendu** :

```
J_s = q/T + ( 3ℓ² / λT₀² ) · q · ∂q/∂z
```

Sans ce terme supplémentaire, la production conserve un terme croisé de signe indéfini et n'est pas
une somme de carrés.

**Piège de dérivation.** Le flux d'entropie doit être pris en `q/T`, non en `q/T₀`. La différence est
exactement le terme qui compense la contribution du gradient de température ; avec `q/T₀` il subsiste
un résidu en `q·∂T/∂z` qui ne s'annule pas.

### 10.4 Système caractéristique

Partie principale du système de Cattaneo en variables `(u, q)` :

```
A = [[0, 1], [a/τ_R, 0]]          valeurs propres  ±√(a/τ_R)
```

Deux valeurs propres réelles distinctes : strictement hyperbolique. Le théorème de Godunov et Mock
relie cette propriété à l'existence d'une entropie convexe.

Implanté : `src/admissibility.py`. Vérifications : `tests/test_admissibility.py`, dont la dérivation
symbolique des deux productions d'entropie.

## 11. Données matériau et limites du cadre de Callaway

Sources dépouillées : Ma, Li, Luo, *Phys. Rev. B* **90**, 035203 (2014) ; Cheng *et al.*,
*Phys. Rev. Materials* **4**, 044602 (2020).

### 11.1 Conductivité mesurée et calculée de l'AlN

| Échantillon | κ à 300 K | Source |
|---|---|---|
| Film MOCVD **sur saphir**, 18 et 22,5 µm | **321 W·m⁻¹·K⁻¹** | Cheng *et al.*, mesure |
| Cristal parfait | 318 W·m⁻¹·K⁻¹ | DFT |
| Substrat massif PVT, Samp_3 | 278 W·m⁻¹·K⁻¹ | Cheng *et al.*, mesure |
| Substrat massif PVT, Samp_4 | 216 W·m⁻¹·K⁻¹ | Cheng *et al.*, mesure |
| Calcul ab initio, dans le plan | 300 W·m⁻¹·K⁻¹ | Ma *et al.* |
| Calcul ab initio, **hors plan** | 286 W·m⁻¹·K⁻¹ | Ma *et al.* |

**Anisotropie réelle : 5 %.** Les films MOCVD sur saphir dépassent la valeur de 285 longtemps
admise comme référence pour l'AlN massif.

Autres valeurs établies : structure wurtzite, quatre atomes par maille, trois branches acoustiques et
neuf optiques. Densité de dislocations des films MOCVD sur saphir, 1,6 × 10⁸ cm⁻². Lacunes
d'aluminium négligeables dans les films MOCVD, 3 × 10¹⁹ à 1,5 × 10²⁰ cm⁻³ dans les substrats PVT.
Mesures de 80 à 480 K.

### 11.2 Correction du contraste d'effusivité

L'hypothèse de travail initiale, une conductivité de film de 60 W·m⁻¹·K⁻¹, était **très pessimiste**.

| Cas | κ | b film | b₃₂ | Γ |
|---|---|---|---|---|
| Valeurs mesurées, film sur saphir | 321 | 27 814 | 0,370 | **+0,46** |
| Substrat PVT le moins bon | 216 | 22 816 | 0,451 | +0,38 |
| Hypothèse initiale | 60 | 12 025 | 0,856 | +0,08 |
| Point aveugle | 44 | 10 298 | 1,000 | 0 |

**Le contraste réel est de l'ordre de 0,46, non de 0,08.** C'est le meilleur régime accessible sur la
carte de conception : l'incertitude y tombe à environ 0,2 %.

**La cécité d'interface n'est donc pas un risque pratique pour l'AlN sur saphir.** Le point aveugle à
44 W·m⁻¹·K⁻¹ exigerait un film très dégradé.

Nuance : les films mesurés font 18 à 22,5 µm. Le libre parcours moyen valant `Λ = 3κ/(Cv) ≈ 67 nm` à
300 K, un film de 500 nm donne `d/Λ ≈ 7,5` : la diffusion aux frontières réduira la conductivité d'un
facteur, non d'un ordre de grandeur.

### 11.3 Lois d'échelle en fréquence, calculées ab initio pour l'AlN wurtzite

```
1/τ_U  ∝  ω³        branches TA et LA
1/τ_N  ∝  ω         branche TA, indépendant de la symétrie du réseau
1/τ_N  ∝  ω²        branche LA, conforme à la prédiction de Herring pour le réseau hexagonal
```

Le rapport `τ_R/τ_N` varie donc comme `ω⁻²` pour les modes transverses. La coordonnée `x` de la carte
des régimes est **fortement dépendante de la fréquence** : un matériau occupe une plage, non un point.

### 11.4 Portée de la diagonale aveugle

Les temps `τ_R` et `τ_ℓ` de Guyer–Krumhansl sont des **moyennes pondérées sur le spectre**. La
condition `τ_R = 1,8 τ_N` porte sur ces moyennes, non sur chaque mode. L'énoncé reste bien défini au
niveau du modèle effectif, mais cette réserve doit accompagner toute confrontation à un matériau réel.

### 11.5 Le cadre de Callaway échoue précisément hors plan

Écart à la solution exacte de l'équation de Boltzmann, AlN à 300 K :

| Direction | RTA | Callaway | Allen modifié |
|---|---|---|---|
| Dans le plan | −11 % | ~0 % au-dessus de 150 K | +10 % |
| **Hors plan** | **−12,3 %** | **−11,6 %** | −8 % |

**Le modèle de Callaway n'apporte qu'une correction de 0,7 % au-dessus de la RTA hors plan**, alors
que la correction nécessaire est de 12,3 %. Dans la direction qui est celle du présent travail, il ne
corrige donc pratiquement rien.

Anisotropie prédite : RTA 7 %, Callaway 19 %, Allen 29 %, contre 5 % en réalité. Le modèle de
Callaway dégrade l'anisotropie au lieu de l'améliorer.

**Raison identifiée par les auteurs.** Entre 50 et 80 THz, de nombreux modes présentent un produit
`v·q` négatif dans la direction hors plan : la vitesse de groupe y est opposée au vecteur d'onde. Le
terme correctif de Callaway change alors de signe et s'annule en moyenne.

### 11.6 Position de l'AlN sur la carte des régimes

Temps de relaxation relevés sur la figure 7 de Ma *et al.*, AlN wurtzite à 300 K, entre 2 et 10 THz.
**Lecture à l'œil sur échelle logarithmique : précision d'un facteur deux.** Elle se valide
toutefois elle-même — les exposants ajustés de `τ_N` sortent à `−2,00` pour LA et `−1,00` pour TA,
exactement les valeurs publiées.

| Branche | ω | τ_U | τ_N | x = τ_U/τ_N |
|---|---|---|---|---|
| LA | 2 THz | 7 × 10⁵ ps | 1,0 × 10⁴ ps | 70 |
| LA | 10 THz | 1,0 × 10⁴ ps | 4 × 10² ps | 25 |
| TA | 2 THz | 4 × 10⁵ ps | 2,0 × 10³ ps | 200 |
| TA | 10 THz | 4 × 10³ ps | 4 × 10² ps | 10 |

**Abscisse de l'AlN à 300 K : `x = 10` à `200`.** La droite de cécité étant à 1,8, le facteur de
sécurité est de 6 au pire. **L'AlN n'est pas près de la cécité à température ambiante.**

La plage d'un facteur 20 sur une seule décade de fréquence confirme quantitativement que la
coordonnée `x` n'est pas une constante du matériau.

**Ordonnée.** Le libre parcours moyen des processus normaux vaut `Λ_N = v·τ_N`, soit **2,4 à 60 µm**
avec `v = 6000 m/s`. Très grand devant un film mince.

| Épaisseur | y = τ_B/τ_N | Position |
|---|---|---|
| 500 nm | 0,008 à 0,21 | très en dessous de la fenêtre |
| 5 µm | 0,08 à 2,1 | à la frontière |
| 20 µm | 0,33 à 8,3 | à la frontière |
| 500 µm | 8,3 à 208 | dans la fenêtre |

**Conséquence : un film d'AlN submicronique est hors de la fenêtre hydrodynamique**, de deux décades.
Les processus normaux y sont trop rares devant la diffusion aux frontières pour que le gaz de phonons
s'équilibre intérieurement. Le terme non local de Guyer–Krumhansl n'y a donc pas de fondement.

### 11.7 Deux libres parcours moyens à ne pas confondre

L'ordonnée de la carte emploie le libre parcours des **processus normaux seuls**. Elle ne dit rien du
caractère balistique du transport, qui dépend du libre parcours **total**.

| Longueur | Expression | Valeur AlN à 300 K | Ce qu'elle décide |
|---|---|---|---|
| Normal | `Λ_N = v τ_N` | 2,4 à 60 µm | possibilité de l'hydrodynamique |
| **Total** | `Λ = 3λ/(ρc v)` | **67 nm** | caractère balistique du transport |

Elles diffèrent de trois ordres de grandeur. Confondre les deux conduit à qualifier de balistique un
film qui ne l'est pas.

**Nombre de Knudsen**, `Kn = Λ/d`, avec `Λ = 67 nm` :

| Épaisseur | Kn | Régime de transport |
|---|---|---|
| 50 nm | 1,33 | balistique |
| 200 nm | 0,33 | transitionnel |
| **500 nm** | **0,13** | **transitionnel** |
| 2 µm | 0,033 | diffusif |

**Un film d'AlN de 500 nm n'est pas balistique : il est transitionnel.** Le libre parcours vaut un
huitième de l'épaisseur. Les phonons y diffusent plusieurs fois en traversant, mais la diffusion aux
frontières reste sensible.

**Conséquence pour le livrable.** Le modèle de diffusion reste applicable, mais **la conductivité
extraite est une valeur apparente, réduite par la diffusion aux frontières et dépendante de
l'épaisseur.** Ce n'est pas une propriété intrinsèque du matériau, et la comparer à une valeur massive
n'a pas de sens. C'est précisément ce que Hoque *et al.* mesurent en faisant varier l'épaisseur de
1,6 à 2440 nm.

**L'épaisseur des films est donc la première caractéristique à établir**, avant la bande de
fréquences : elle décide de ce que le livrable peut annoncer.

Implanté : `Sample.mean_free_path`, `Sample.knudsen_number`, `Sample.transport_regime`,
`regime_report()`.

Figure : `figures/05_regime_map.png`.

### 11.8 Ce que les sources ne fournissent pas

Le modèle de Cheng *et al.* **n'inclut pas les processus normaux** :

```
1/τ_C = 1/τ_U + 1/τ_M + 1/τ_B
1/τ_U = B T ω² exp(−C/T)
1/τ_M = ( V ω⁴ / 4π v³ ) Σ x_i (ΔM_i/M)²
1/τ_B = v / d
```

Il s'agit d'une RTA avec règle de Matthiessen. Il ne peut donc **pas fournir `τ_N`**, donc pas
l'abscisse de la carte des régimes. Les auteurs attribuent eux-mêmes leurs écarts à basse température
aux limites du modèle de Callaway.

Pour obtenir `τ_N`, il faut soit les données de la figure 7 de Ma *et al.*, soit le modèle de
Debye–Callaway modifié de Morelli, qui exprime les coefficients par les paramètres de Grüneisen.

## 12. Contraintes instrumentales établies

### 12.1 Fréquence caractéristique du film

```
f_c = 1 / (2π ξ₁²)        ξ₁ = e/√a
```

En dessous, la réponse est dominée par le substrat et le film n'est pas accessible. Le confinement de
l'onde dans le film exige environ deux décades au-dessus.

**Film d'AlN de 500 nm : `f_c ≈ 16 MHz`.** Un banc travaillant en kilohertz ne voit pas le film.

Implanté : `Sample.characteristic_frequency`.

### 12.2 Effusivité aveugle

```
Γ = (1 − b₃₂)/(1 + b₃₂)        b₃₂ = b_substrat / b_film
```

`Γ` est le coefficient de réflexion de l'onde thermique à l'interface. Il s'annule quand les deux
effusivités sont égales, et **l'interface devient alors strictement invisible** : la réponse de face
avant est celle d'un milieu semi-infini, pour toute épaisseur et tout contraste de diffusivité.
Vérifié numériquement à 10⁻¹² près.

**AlN sur saphir : l'interface est aveugle pour une conductivité de film de 44,0 W·m⁻¹·K⁻¹.**
L'incertitude passe de 8 % à 2300 % pour un écart de 1 W·m⁻¹·K⁻¹ autour de cette valeur.

**Mais ce point est loin des valeurs réelles.** Les films d'AlN sur saphir mesurent 321 W·m⁻¹·K⁻¹, et
les substrats les plus dégradés 216, soit `Γ` entre 0,38 et 0,46. Voir section 11.2. La cécité
d'interface est donc un cas limite du formalisme, non un risque pratique pour ce système. Elle le
redeviendrait pour un film fortement dégradé, ou pour un autre couple film-substrat.

Implanté : `Sample.reflection_coefficient`, `Sample.blind_film_effusivity`, `contrast_report()`.

### 12.3 Seuil de mesurabilité du temps de relaxation

```
ω_max · τ ≥ 1        soit        τ ≥ 1/(2π f_max)
```

En dessous, l'incertitude relative croît comme l'inverse de `ω_max·τ`. La sensibilité
`dφ/d ln τ = ½ ωτ/(1+(ωτ)²)` est maximale exactement en `ωτ = 1`, où elle vaut 1/4 radian, et décroît
symétriquement de part et d'autre en échelle logarithmique.

**Plafond d'un banc FDTR à 200 MHz : `τ_min = 8,0 × 10⁻¹⁰ s`.** Cette valeur ne dépend que de
l'instrument, pas du matériau.

**Conséquence pour l'AlN :** les temps de relaxation phononiques y sont de l'ordre de la picoseconde
à quelques dizaines de picosecondes, soit deux à trois ordres de grandeur en dessous du plafond. **Le
temps de relaxation de l'AlN n'est pas mesurable par thermoréflectance en domaine fréquentiel.** Ce
constat s'appuie sur la bande que Camacho de la Rosa et al. (2025) annoncent eux-mêmes.

Le même seuil apparaît par trois voies indépendantes : l'argument de l'énergie de l'analogie
spectrale vaut −45° en `ωτ = 1` ; la phase mesurable y franchit sa médiane à −22,5° ; et le plafond
instrumental impose `ω_max·τ ≥ 1`. Figure `figures/02_omega_tau_threshold.png`.

### 12.4 Diagonale aveugle des deux temps

```
τ_R = τ_ℓ        soit, en grandeurs microscopiques,        τ_R = 1,8 τ_N
```

Quand les deux temps coïncident, les racines des deux combinaisons se simplifient et la réponse est
**exactement celle d'un milieu de Fourier**, quelle que soit leur valeur commune. Vérifié à 10⁻¹³
sur quatre décades de temps et huit décades de fréquence.

Un tel matériau est thermiquement indiscernable d'un matériau de Fourier. Aucune mesure thermique ne
peut le distinguer.

**Critère de séparabilité**, au bruit de référence : les deux temps doivent différer d'au moins
**16 %** pour une incertitude inférieure à 10 %, et d'au moins **25 %** pour une incertitude
inférieure à 1 %.

Cette cécité diffère de celle de la section 10.2 par sa nature. Celle de l'interface dépend du choix
de substrat, donc elle se contourne. Celle-ci est **intrinsèque au matériau** : aucun changement de
dispositif ne la lève. Elle dépend en revanche de la température, `τ_N` et `τ_R` suivant des lois
différentes — un cristal peut donc traverser la condition en refroidissant.

Figure `figures/04_guyer_krumhansl_blindness.png`.

### 12.5 Débordement numérique

```
|√P · ξ| < 700
```

Au-delà, `cosh` et `sinh` débordent en double précision. Sous Cattaneo, `P ≈ τω²` croît
quadratiquement, donc l'argument croît **linéairement** en fréquence au lieu de la racine : la bande
exploitable se rétrécit comme l'inverse du carré de la fréquence maximale quand `τ` augmente.

Limite de la représentation matricielle, non de la physique. Une formulation à facteur exponentiel
extrait la lèverait.

## 13. Inversion numérique de Laplace

**Méthode initiale : Gaver–Stehfest**, Maillet appendice 1.1 page 28. **Repli : De Hoog**, employée
par Krapez.

Erreur relative mesurée à `n = 12`, selon la décroissance du signal depuis son maximum :

| Signal restant | Erreur relative |
|---|---|
| au-dessus de 10⁻¹ | 5 × 10⁻⁴ |
| au-dessus de 10⁻² | 7 × 10⁻³ |
| au-dessus de 10⁻³ | 10⁻¹ |

**Règle : l'inversion tient sur environ deux décades de décroissance.** Un ajustement sur une queue
au-delà porte sur du bruit numérique. La précision n'est pas monotone en `n` : elle s'améliore
jusqu'à 12 puis se dégrade.

Validation obligatoire avant tout usage sur bicouche : comparaison à la solution analytique du mur
homogène.

## 14. Validation sur données synthétiques

Valider une inversion avec le modèle même qui a engendré les données porte un nom : **crime inverse**
(Krapez 2023, section IV).

Ce n'est pas illégitime — c'est le seul moyen de vérifier qu'un estimateur retrouve une vérité
connue — mais il faut le nommer et dire ce qu'il ne prouve pas. Il ne teste ni l'adéquation du modèle
à la réalité, ni la robustesse aux effets absents du modèle : rugosité, conductance d'interface,
phase de référence du montage, anisotropie.

Une barre d'erreur n'est crédible que si elle a été confrontée à la dispersion observée sur de
nombreuses réalisations de bruit. Voir `test_fisher_uncertainties_match_monte_carlo`.

### Validation externe

**Toute loi physique implémentée exige au moins une comparaison à une solution de référence prise
hors du code.**

Les tests de cohérence interne — déterminant unimodulaire, limite homogène, composition, découpage —
ne détectent pas une erreur cohérente avec elle-même. L'épisode consigné en section 7.4 en donne un
cas : une erreur de coefficient de flux a survécu à cent treize tests, et seule la confrontation à
la forme fermée d'un article publié l'a révélée.

Références de validation en usage :

| Loi | Référence externe |
|---|---|
| Fourier, mur homogène | Maillet et al. §1.3.4, et solution analytique `℘/(b√p)` |
| Fourier, impulsionnel | `Q/(b√(πt))` |
| Cattaneo, semi-infini | Camacho de la Rosa et al. (2025), phase `−45° + ½·arctan(ωτ)` |
| Guyer–Krumhansl | réduction exacte à Cattaneo et à Fourier |

## 15. Tests unitaires associés

| Test | Critère |
|---|---|
| Mur homogène | une expression graduée à effusivités égales redonne l'expression de la section 2 |
| Déterminant | `AD − BC = 1`, nécessaire mais non suffisant |
| Comparaison terme à terme | chaque entrée comparée séparément, pour détecter un écart de convention |
| Découpage | une couche coupée en deux redonne la couche entière |
| Fermeture semi-infinie | film et substrat de mêmes propriétés donnent `℘/(b√p)` exactement |
| Interface aveugle | effusivités égales, la réponse est celle du semi-infini pour toute épaisseur |
| Inversion de Laplace | écart borné à la solution analytique sur la plage utile |
| Incertitudes | accord entre prédiction de Fisher et dispersion Monte-Carlo |
| Relation d'Euler | résidu nul sous Fourier, sous Cattaneo et sous Guyer–Krumhansl |
| Terme en `ℓ` requis | omettre la dérivée en `ℓ` brise l'identité d'un facteur mille |
| Limite Fourier | `τ_R = τ_ℓ = 0` redonne le modèle de Fourier à la précision machine |
| Limite Cattaneo | `τ_ℓ = 0` redonne la forme de Cattaneo |
| Diagonale aveugle | `τ_R = τ_ℓ` redonne la réponse de Fourier exactement |
| Validation externe | phase comparée à la forme fermée de Camacho de la Rosa et al. |
| Dispersion | exposant de `k` en `ω` : 0,5 pour Fourier, 1,0 pour Cattaneo, 0,5 pour Guyer–Krumhansl |
| Vitesse de Cattaneo | sature à `√(a/τ_R)`, soit `v/√3` en modèle de Debye |
| Production d'entropie | somme de carrés, dérivée symboliquement depuis les bilans |
| Flux d'entropie étendu | coefficient `3ℓ²/(λT₀²)`, obtenu en annulant le terme croisé |
| Refus explicite | couche graduée avec relaxation : `NotImplementedError`, pas un résultat faux |
