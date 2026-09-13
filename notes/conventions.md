# Conventions retenues

Ce fichier fixe les conventions du projet. Il fait autorité sur toute expression employée dans `src/`.
Toute source extérieure doit être convertie vers ces conventions avant usage.

Dernière mise à jour : 13 septembre 2026.

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

## 10. Contraintes instrumentales établies

### 10.1 Fréquence caractéristique du film

```
f_c = 1 / (2π ξ₁²)        ξ₁ = e/√a
```

En dessous, la réponse est dominée par le substrat et le film n'est pas accessible. Le confinement de
l'onde dans le film exige environ deux décades au-dessus.

**Film d'AlN de 500 nm : `f_c ≈ 16 MHz`.** Un banc travaillant en kilohertz ne voit pas le film.

Implanté : `Sample.characteristic_frequency`.

### 10.2 Effusivité aveugle

```
Γ = (1 − b₃₂)/(1 + b₃₂)        b₃₂ = b_substrat / b_film
```

`Γ` est le coefficient de réflexion de l'onde thermique à l'interface. Il s'annule quand les deux
effusivités sont égales, et **l'interface devient alors strictement invisible** : la réponse de face
avant est celle d'un milieu semi-infini, pour toute épaisseur et tout contraste de diffusivité.
Vérifié numériquement à 10⁻¹² près.

**AlN sur saphir : l'interface est aveugle pour une conductivité de film de 44,0 W·m⁻¹·K⁻¹**, valeur
qui tombe au milieu de la plage plausible d'un film mince. L'incertitude passe de 8 % à 2300 % pour
un écart de 1 W·m⁻¹·K⁻¹.

Implanté : `Sample.reflection_coefficient`, `Sample.blind_film_effusivity`, `contrast_report()`.

### 10.3 Seuil de mesurabilité du temps de relaxation

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

### 10.4 Diagonale aveugle des deux temps

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

### 10.5 Débordement numérique

```
|√P · ξ| < 700
```

Au-delà, `cosh` et `sinh` débordent en double précision. Sous Cattaneo, `P ≈ τω²` croît
quadratiquement, donc l'argument croît **linéairement** en fréquence au lieu de la racine : la bande
exploitable se rétrécit comme l'inverse du carré de la fréquence maximale quand `τ` augmente.

Limite de la représentation matricielle, non de la physique. Une formulation à facteur exponentiel
extrait la lèverait.

## 11. Inversion numérique de Laplace

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

## 12. Validation sur données synthétiques

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

## 13. Tests unitaires associés

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
| Refus explicite | couche graduée avec relaxation : `NotImplementedError`, pas un résultat faux |
