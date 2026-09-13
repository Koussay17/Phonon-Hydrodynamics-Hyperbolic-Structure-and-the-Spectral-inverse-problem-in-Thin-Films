# Conventions retenues

Ce fichier fixe les conventions du projet. Il fait autorité sur toute expression employée dans `src/`.
Toute source extérieure doit être convertie vers ces conventions avant usage.

Dernière mise à jour : 12 septembre 2026.

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

## 7. Loi constitutive non-Fourier

**Substitution retenue, loi de Cattaneo à un seul temps de relaxation :**

```
p  →  P = p (1 + τ p)
```

En régime modulé, `P = iω − τω²` : une partie réelle apparaît, croissant comme le carré de la
fréquence. Le groupement sans dimension est `ωτ`.

**Rien d'autre ne change.** La coordonnée `ξ`, l'effusivité `b` et le potentiel `V` gardent leurs
définitions de Fourier. Le quadripôle est celui de Fourier évalué en `P`.

Vérification croisée : Camacho de la Rosa et al. (2025) obtiennent la même substitution par une voie
indépendante — fonction de Green tridimensionnelle et transformée de Hankel — sous la forme
`σ²_cv = σ²_f (1 + iωτ)`.

Chaque couche est évaluée à son propre `P` : un temps de relaxation dans le film n'affecte pas le
substrat et réciproquement.

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

### 9.2 Sous Cattaneo

**L'invariance survit intacte.** `τ` est un temps, non une longueur : un changement d'échelle en
profondeur ne l'affecte pas, donc `P` est invariant, donc le quadripôle l'est aussi.

Vérifié numériquement : le résidu d'Euler reste au niveau de l'erreur de troncature quel que soit `τ`,
jusqu'à `ωτ` de plusieurs centaines. Voir `tests/test_cattaneo.py`.

`τ` entre en revanche comme un paramètre authentiquement indépendant : son estimation ne dégrade ni
celle de `b` ni celle de `ξ₁`.

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

### 10.4 Débordement numérique

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
| Relation d'Euler | résidu nul sous Fourier, et sous Cattaneo |
| Limite Fourier | `τ = 0` redonne le modèle de Fourier à la précision machine |
