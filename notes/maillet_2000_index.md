# Maillet et al. (2000) — index de consultation

**Référence.** D. Maillet, S. André, J.-C. Batsale, A. Degiovanni, C. Moyne,
*Thermal Quadrupoles: Solving the Heat Equation through Integral Transforms*,
John Wiley & Sons, New York, 2000.

Ouvrage de référence du formalisme quadripolaire. Il se consulte par entrée, il ne se lit pas
linéairement. La pagination ci-dessous correspond à l'exemplaire utilisé et doit être vérifiée en cas
d'édition différente.

Cet index est destiné à être complété au fil du travail.

---

## 1. Quadripôle du mur homogène

| Section | Page | Contenu |
|---|---|---|
| §3.3.1 | 73 | *The Notion of a Thermal Quadrupole: the Passive 'Wall'* — dérivation de référence |
| §1.3.4 | 10 | *Case d: A One-layer Slab in Transient Transfer* — forme courte, directement utilisable |
| §4.3 | 159 | Extension au traitement bidimensionnel et tridimensionnel |

Expression donnée en §1.3.4, avec `k = √(p/a)`, `e` l'épaisseur et `S` la section :

```
A = D = cosh(ke)
B = sinh(ke) / (λ k S)
C = λ k S sinh(ke)
```

Correspondance avec la coordonnée transformée : `ke = √p · ξ₁` et `λk = b√p`, d'où l'identité des
deux écritures au facteur `S` près. Voir `conventions.md`.

## 2. Conventions de signe et vecteur d'état

| Section | Page | Contenu |
|---|---|---|
| §1.3.1 | 4 | Sens entrée vers sortie du vecteur d'état, figures 1.1 à 1.3 |
| §2.2.1 à §2.2.4 | 37–43 | Flux, densité de flux, loi de Fourier ; signe et orientation de l'axe |
| §2.3.2 | 44 | Conditions aux limites |
| §2.3.3 | 51 | Milieux multiples et conditions d'interface ; orientation aux raccords |
| §2.5 | 61 | Classification des types de problèmes ; mise en forme du vecteur `[θ φ]ᵀ` |

Les sections 2.3.3 et 2.5 sont celles dont dépend la validité du chaînage matriciel.

## 3. Fermeture par impédance sur substrat semi-infini

Aucune entrée de la table des matières ne porte ce nom. Les points d'entrée effectifs sont :

| Section | Page | Contenu |
|---|---|---|
| §1.4.2, remarque 2 | 15 | Passage à la limite d'épaisseur infinie donnant `Z = 1/(b√p)` par unité de surface. **Formule de fermeture elle-même.** |
| §3.3.3 | 84 | Quadripôle associé à une condition de Fourier ; forme générique de la matrice de fermeture |
| §3.5.1 | 102 | Terminaison par impédance sur un élément fluide ; transposable |
| §5.4.2 | 192 | Plaque semi-infinie ou finie en régime périodique établi |
| §6.1 | 211–240 | Constriction des lignes de flux ; impédance de constriction si la section de contact est réduite |

La valeur de la remarque 2 coïncide avec l'impédance employée par Krapez à son équation (29).

## 4. Inversion numérique de Laplace

| Section | Page | Contenu |
|---|---|---|
| Chapitre 9 | 333 | Ensemble du chapitre |
| §9.2.4 | 337 | Table 9.1 des transformées usuelles et inversions explicites |
| §9.3 | 340 | Transformée inverse |
| §9.3.2 | 342 | **Méthode de Gaver–Stehfest** |
| §9.4 | 345 | Comparaison des méthodes |
| §9.5 | 351 | Recommandations pratiques |
| Appendice 1.1 | 28 | Coefficients de Stehfest et implémentation MATLAB |

Krapez emploie pour sa part la méthode de De Hoog. Stehfest est plus simple à implémenter mais se
dégrade sur les réponses oscillantes et requiert une arithmétique en précision étendue. Toute
implémentation doit être validée contre le cas analytique du mur homogène avant application à un
bicouche.

---

## À compléter

- Résistance thermique de contact et quadripôle d'interface associé
- Sources internes, localisées et distribuées
- Régime périodique établi : formulation complète
- Estimation de paramètres et sensibilités, si le chapitre correspondant existe
