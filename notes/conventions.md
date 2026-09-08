# Conventions retenues

Ce fichier fixe les conventions du projet. Il fait autorité sur toute expression employée dans `src/`.
Toute source extérieure doit être convertie vers ces conventions avant usage.

Dernière vérification : 8 septembre 2026.

---

## 1. Grandeur de flux

**Convention retenue : densité de flux, par unité de surface, en W·m⁻².**

Deux formalismes coexistent dans la littérature consultée.

| Source | Grandeur | Unité | Section apparente |
|---|---|---|---|
| Maillet et al. (2000) | flux total | W | oui, facteur `S` |
| Krapez (2018) | densité de flux | W·m⁻² | non |
| Présent travail | densité de flux | W·m⁻² | non |

Motifs du choix : le problème traité est unidimensionnel, la section n'y joue aucun rôle physique ;
et les comparaisons avec Krapez se font alors sans conversion.

**Règle d'usage.** Toute expression tirée de Maillet portant un facteur `S` s'emploie avec `S = 1`.

## 2. Vérification effectuée sur le mur homogène

Maillet, §1.3.4 page 10, avec `k = √(p/a)`, `e` l'épaisseur, `S` la section :

```
A = D = cosh(ke)
B = sinh(ke) / (λ k S)
C = λ k S sinh(ke)
```

Expression obtenue par la voie transformée, avec `ξ₁ = e/√a` :

```
A = D = cosh(√p ξ₁)
B = sinh(√p ξ₁) / (b √p)
C = b √p sinh(√p ξ₁)
```

Correspondance des notations, par les identités `λ/√a = b` et `ξ₁ = e/√a` :

```
k e  = √(p/a) · e = √p · e/√a = √p · ξ₁
λ k  = (λ/√a) · √p = b √p
```

**Conclusion.** Les deux expressions sont identiques au seul facteur `S` près, qui figure au
dénominateur de `B` et au numérateur de `C`. Aucun autre écart de convention n'a été relevé :
orientation de l'axe, sens du vecteur d'état et signe du flux coïncident.

## 3. Point de vigilance sur le déterminant

Le déterminant vaut 1 dans les deux formalismes :

```
A D − B C = cosh² − [sinh/(λkS)] · [λkS sinh] = cosh² − sinh² = 1
```

Le facteur `S` se simplifie. **Le test du déterminant ne détecte donc pas un mélange de conventions.**
Un produit matriciel associant une matrice avec `S` et une matrice sans aboutit, satisfait le test, et
donne un résultat faux d'un facteur `S`.

Ce test reste nécessaire mais n'est pas suffisant. La validation d'une expression quadripolaire exige
en outre la comparaison terme à terme avec le cas homogène.

## 4. Impédance de fermeture

**Substrat semi-infini, par unité de surface :**

```
Z = 1 / (b √p)
```

Vérification croisée : Maillet §1.4.2 remarque 2 page 15, et Krapez équation (29). Les deux sources
donnent la même expression, ce qui confirme que le seul écart entre les formalismes portait sur `S`.

Réponse en surface d'un système revêtement sur substrat, avec `h` le coefficient d'échange en face
avant et `℘` la densité de puissance de la source :

```
θ₀ = ℘ (A Z + B) / (C Z + D + h (A Z + B))
```

Le cas adiabatique correspond à `h = 0`.

## 5. Coordonnée et grandeurs dérivées

```
ξ(z) = ∫₀ᶻ du / √(a(u))        [ξ] = s^(1/2)
a = λ / ρc                      [a] = m²·s⁻¹
b = √(λ ρc)                     [b] = W·s^(1/2)·m⁻²·K⁻¹
λ = b √a          ρc = b / √a
```

Forme normale et potentiel :

```
d²ψ/dξ² − (V(ξ) + p) ψ = 0
V = s″ / s        s = b^(±1/2)
```

Exposant positif pour la formulation en température, négatif pour la formulation en flux.

## 6. Inversion numérique de Laplace

**Méthode initiale : Gaver–Stehfest**, coefficients et implémentation dans Maillet, appendice 1.1
page 28.

**Méthode de repli : De Hoog**, employée par Krapez.

Critère de bascule : Stehfest se dégrade sur les réponses oscillantes et requiert une arithmétique en
précision étendue. Si l'écart au cas analytique du mur homogène croît aux temps courts, passer à
De Hoog.

**Validation obligatoire avant tout usage sur bicouche :** comparaison à la solution analytique du mur
homogène.

## 7. Tests unitaires associés

| Test | Critère |
|---|---|
| Mur homogène | Une expression graduée avec effusivités égales aux deux faces redonne l'expression de la section 2 |
| Déterminant | `A D − B C = 1`, nécessaire mais non suffisant |
| Comparaison terme à terme | Chaque entrée comparée séparément au cas homogène, pour détecter un écart de convention |
| Inversion de Laplace | Écart à la solution analytique du mur homogène borné sur toute la plage temporelle utile |
