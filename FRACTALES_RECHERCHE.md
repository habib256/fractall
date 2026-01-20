# Recherche Bibliographique : Nouveaux Types de Fractales

Ce document consolide les résultats de la Phase 1 de recherche bibliographique pour les nouveaux types de fractales à implémenter dans fractall.

**Date de recherche** : 2024
**Statut** : Phase 1 complétée

---

## 1. Perpendicular Burning Ship

### Formule
```
z₀ = 0
z_{n+1} = (Re(z) - i·|Im(z)|)² + c
```

En termes de parties réelle et imaginaire :
```
x_{n+1} = x_n² - |y_n|² + c_x
y_{n+1} = -2·x_n·|y_n| + c_y
```

### Paramètres par défaut recommandés
- **Domaine** : `[-2.5, 1.5] × [-1.5, 1.5]` (similaire à Mandelbrot)
- **iterationMax** : 1000-5000 (selon niveau de détail souhaité)
- **bailout** : 4.0 (seuil d'échappement standard)
- **seed** : `(0, 0)` (Mandelbrot-style)
- **zoomfactor** : 8 (comme Burning Ship)

### Références
- **Source originale** : Kalle's Fractaler 2 (kf2)
- **Code de référence** : 
  - GitHub smurfix/kf2 : Formule "Perpendicular Burning Ship"
  - FractalShades (gbillotey/Fractalshades) : Implémentation Python
- **Documentation** : 
  - GitHub kf2 formula listing
  - FractalShades documentation

### Notes techniques
- **GMP** : Oui (compatible avec précision arbitraire)
- **SIMD** : Oui (opérations simples sur doubles)
- **Particularités** : 
  - Seule la partie imaginaire est pliée (abs), contrairement au Burning Ship standard
  - Le signe négatif dans y_{n+1} crée l'effet "perpendicular"
  - Compatible avec l'architecture escape-time existante
  - Structure similaire à Burning Ship (type 13), donc facile à intégrer

### Images de référence
- Formes géométriques nettes avec structures symétriques
- Zones intéressantes : régions centrales avec motifs répétitifs

---

## 2. Celtic Fractal

### Formule
```
z₀ = 0
z_{n+1} = |Re(z²)| + i·Im(z²) + c
```

En termes de parties réelle et imaginaire :
```
u = x² - y²        // partie réelle de z²
v = 2·x·y          // partie imaginaire de z²
x_{n+1} = |u| + c_x
y_{n+1} = v + c_y
```

### Paramètres par défaut recommandés
- **Domaine** : `[-2.0, 1.0] × [-1.5, 1.5]` (similaire à Mandelbrot)
- **iterationMax** : 1000-5000
- **bailout** : 4.0
- **seed** : `(0, 0)` (Mandelbrot-style)
- **zoomfactor** : 8

### Variantes
- **Wide Celtic** : `z_{n+1} = |Re(z·|Re(z)|) - Im(z)²| + i·2·Re(z)·Im(z) + c`
- **Celtic Buffalo** : Application de abs() aux deux parties

### Références
- **Source originale** : Communauté fractale (UltraFractal, DeviantArt)
- **Code de référence** :
  - Reddit r/fractals : Discussions sur formules Celtic
  - UltraFractal : Formules documentées
- **Documentation** :
  - Reddit discussions avec formules précises
  - UltraFractal formula database

### Notes techniques
- **GMP** : Oui
- **SIMD** : Oui
- **Particularités** :
  - Modification simple : abs() appliqué uniquement à la partie réelle après le carré
  - Crée des motifs entrelacés caractéristiques
  - Très populaire dans la communauté fractale
  - Implémentation directe compatible avec escape-time

### Images de référence
- Motifs celtiques entrelacés
- Symétrie complexe avec structures répétitives

---

## 3. Alpha Mandelbrot (Nested/Composite)

### Formule
```
z₀ = 0
z_{n+1} = z² + (z² + c)² + c
```

En termes de parties réelle et imaginaire :
```
m = z² + c
z_{n+1} = z² + m² + c
```

### Paramètres par défaut recommandés
- **Domaine** : `[-2.5, 1.5] × [-1.5, 1.5]` (similaire à Mandelbrot)
- **iterationMax** : 500-2000 (divergence plus rapide)
- **bailout** : 4.0 ou plus élevé (peut diverger rapidement)
- **seed** : `(0, 0)` (Mandelbrot-style)
- **zoomfactor** : 8

### Références
- **Source originale** : Reddit r/fractals (discussion 2023)
- **Code de référence** :
  - Reddit : Post original avec formule exacte
  - Communauté fractale : Variantes "alpha fractals"
- **Documentation** :
  - Reddit discussion : "I've created a new kind of fractal variant called Alpha Mandelbrot"

### Notes techniques
- **GMP** : Oui (pour zooms profonds)
- **SIMD** : Oui (opérations standard)
- **Particularités** :
  - Nested/composite : double carré avant ajout de c
  - Divergence plus rapide que Mandelbrot standard
  - Structures multi-couches avec auto-similarité renforcée
  - Peut nécessiter ajustement du bailout selon zoom

### Images de référence
- Structures superposées avec détails supplémentaires
- Bulges et distortions par rapport à Mandelbrot standard

---

## 4. Pickover Stalks / Biomorphs

### Formule de base
```
z₀ = 0 (ou variante)
z_{n+1} = z² + c  (ou autres formules : sin(z), exp(z), etc.)
```

### Système d'orbit trap
Pendant l'itération, calculer à chaque étape :
```
distance_trap = min(|Re(z)|, |Im(z)|)
trap_min = min(trap_min, distance_trap)  // minimum sur toutes les itérations
```

### Coloration
Utiliser `log(trap_min)` ou `trap_min / divisor` pour la coloration.
Les points dont l'orbite s'approche des axes (Re=0 ou Im=0) produisent des "stalks" brillants.

### Paramètres par défaut recommandés
- **Domaine** : `[-2.0, 1.0] × [-1.5, 1.5]`
- **iterationMax** : 100-1000
- **bailout** : 4.0 ou 100.0 (selon variante)
- **seed** : `(0, 0)` pour Mandelbrot-style
- **zoomfactor** : 4-8
- **transformation vector** : `(0, 0)` par défaut (peut être ajusté pour décaler le trap)
- **colour_divisor** : 0.03 (ajuste l'épaisseur des stalks)

### Variantes de formule
- **Mandelbrot** : `z_{n+1} = z² + c`
- **Biomorphs classiques** : `z_{n+1} = sin(z) + exp(z) + c`
- **Autres** : `z_{n+1} = z^z + z^5 + c`, etc.

### Références
- **Source originale** : Clifford A. Pickover (Pickover stalks)
- **Code de référence** :
  - Paul Bourke : Articles et exemples
  - Mitch Reid : Implémentations avec orbit trap
  - Paul Nylander : Formules et paramètres
- **Documentation** :
  - Wikipedia : "Pickover stalk"
  - Paul Bourke's Fractals : Articles sur biomorphs
  - Nylander's blog : Formules détaillées

### Notes techniques
- **GMP** : Oui (si zoom profond)
- **SIMD** : Oui (calculs standard)
- **Particularités** :
  - Nécessite un système d'orbit trap spécialisé
  - Algorithme différent : pas seulement escape-time, mais suivi de distance minimale
  - Coloration basée sur trap_min plutôt que iteration count
  - Compatible avec escape-time mais nécessite modification de la fonction de coloration

### Images de référence
- Formes biomorphiques organiques
- Structures en croix (stalks) le long des axes
- Formes biologiques avec détails fins

---

## 5. Nova Fractal

### Formule
```
z₀ = valeur initiale (souvent point critique ou 0)
z_{n+1} = z - a·(p(z)/p'(z)) + c
```

Pour le polynôme classique `p(z) = z^p - 1` :
```
z_{n+1} = z - a·((z^p - 1)/(p·z^(p-1))) + c
```

Pour p=3 (cubique) :
```
z_{n+1} = z - a·((z³ - 1)/(3·z²)) + c
```

### Paramètres par défaut recommandés
- **Domaine** : `[-3.0, 3.0] × [-2.0, 2.0]` (domaine plus large)
- **iterationMax** : 100-1000
- **bailout** : Variable, souvent 4.0-100.0
- **a (relaxation)** : `(1.0, 0.0)` par défaut (peut être complexe)
- **p (exposant)** : 3 (cubique) ou autre entier
- **seed** : Point critique ou `(0, 0)`
- **zoomfactor** : 4

### Paramètres ajustables
- **a (relaxation)** : Contrôle la vitesse de convergence
  - `a = 1` : Convergence standard
  - `a < 1` : Convergence ralentie
  - `a > 1` : Convergence accélérée, peut créer chaos
- **p (exposant)** : Détermine la symétrie
  - Entier : Symétrie régulière (p-1 lobes)
  - Non-entier : Distorsions et morphing
- **c (seed)** : 
  - Julia mode : c fixe, z₀ varie
  - Mandelbrot mode : c varie, z₀ = point critique

### Références
- **Source originale** : Paul Derbyshire (milieu des années 1990)
- **Code de référence** :
  - UltraFractal : Formule "Nova" standard
  - HPDZ.net : Exemples et galeries
- **Documentation** :
  - UltraFractal help : Documentation complète
  - HPDZ.net : Galeries avec paramètres
  - Wikipedia : "Newton fractal"

### Notes techniques
- **GMP** : Oui (pour précision dans calculs de dérivée)
- **SIMD** : Partiel (calculs complexes avec divisions)
- **Particularités** :
  - Nécessite calcul de dérivée polynomiale
  - Paramètre de relaxation `a` ajoute flexibilité
  - Peut produire des structures très différentes selon paramètres
  - Mode Julia et Mandelbrot possibles
  - Structures en spirales et antennes caractéristiques

### Images de référence
- Bras spiralés élégants
- Motifs répétitifs fins avec structures en antennes
- Zoom profond révèle détails supplémentaires

---

## 6. Multibrot (Puissances Non-entières)

### Formule
```
z₀ = 0
z_{n+1} = z^d + c
```

où `d` est un nombre réel (non-entier), par exemple 2.5, 3.7, etc.

### Calcul de z^d avec branch cuts
```
z = r·e^(iθ)  où r = |z|, θ = arg(z) ∈ (-π, π]
z^d = r^d · e^(i·d·θ)
```

Utiliser la branche principale (principal branch) :
```
log(z) = ln(r) + i·θ  (avec θ ∈ (-π, π])
z^d = exp(d·log(z))
```

### Paramètres par défaut recommandés
- **Domaine** : `[-2.5, 1.5] × [-1.5, 1.5]` (ajuster selon d)
- **iterationMax** : 1000-5000
- **bailout** : 4.0
- **d (exposant)** : 2.5 (exemple de morphing entre 2 et 3)
- **seed** : `(0, 0)`
- **zoomfactor** : 8

### Gestion des branch cuts
- **Branch cut standard** : Le long de l'axe réel négatif
- **Principal branch** : `arg(z) ∈ (-π, π]`
- **Discontinuité** : À la traversée du branch cut, `arg(z)` saute de 2π
- **Multi-branch rendering** : Possible pour visualiser toutes les branches (q branches pour d = p/q)

### Références
- **Source originale** : Généralisation classique de Mandelbrot
- **Code de référence** :
  - Paul Bourke : Articles et exemples
  - Mitch Reid : Animations morphing
  - GitHub : Implémentations diverses
- **Documentation** :
  - Paul Bourke's Fractals : "Mandelbrot Power"
  - Articles sur complex exponentiation et branch cuts

### Notes techniques
- **GMP** : Oui (important pour précision avec puissances non-entières)
- **SIMD** : Partiel (calculs complexes avec exponentielles)
- **Particularités** :
  - Nécessite gestion explicite des branch cuts
  - Calcul de z^d via exponentiation complexe
  - Morphing fluide entre formes entières
  - Symétries variables selon d
  - Discontinuités visuelles possibles au branch cut

### Images de référence
- Transitions fluides entre formes entières
- Symétries variables selon l'exposant
- Détails fins avec structures répétitives

---

## 7. Lambda Fractal (Suite Logistique Complexe)

### Formule
```
z₀ = valeur initiale (souvent dans [0, 1] ou complexe)
z_{n+1} = λ·z·(1 - z)
```

où `λ` est un paramètre complexe.

### Relation avec Mandelbrot
La suite logistique complexe est conjuguée à `z² + c` via transformation :
```
c = λ(2 - λ)/4
```

### Paramètres par défaut recommandés
- **Domaine** : `[0, 4] × [-2, 2]` (domaine du paramètre λ)
- **iterationMax** : 1000-5000
- **bailout** : Variable (peut être 4.0 ou plus)
- **λ (paramètre)** : Variable selon pixel (Mandelbrot-style)
- **z₀** : `0.5` ou `(0.5, 0)` (valeur initiale typique)
- **zoomfactor** : 4-8

### Références
- **Source originale** : Généralisation classique de la suite logistique
- **Code de référence** :
  - Wikipedia : Articles sur logistic map
  - Recherche académique : Articles sur complex logistic map
- **Documentation** :
  - Wikipedia : "Logistic map"
  - Articles scientifiques sur fractales logistiques

### Notes techniques
- **GMP** : Oui (si nécessaire)
- **SIMD** : Oui (opérations simples)
- **Particularités** :
  - Formule simple mais comportement complexe
  - Relation avec Mandelbrot via conjugaison
  - Formes organiques caractéristiques
  - Complémentaire à Lyapunov (type 17)
  - Domaine de paramètres différent de Mandelbrot standard

### Images de référence
- Formes organiques avec structures chaotiques
- Connexions avec théorie du chaos
- Patterns similaires à Lyapunov mais avec formule différente

---

## 8. Anti-Buddhabrot

### Algorithme
Contrairement au Buddhabrot (type 16) qui trace les trajectoires des points qui **s'échappent**, l'Anti-Buddhabrot trace les trajectoires des points qui **ne s'échappent pas** (restent dans le Mandelbrot set).

### Étapes algorithmiques
1. Échantillonnage de points `c` dans le Mandelbrot set (intérieur)
2. Pour chaque `c`, itérer `z_{n+1} = z² + c` jusqu'à détection de cycle périodique
3. Tracer toutes les valeurs `z` visitées (ou seulement le cycle limite)
4. Incrémenter compteur de densité pour chaque pixel traversé
5. Coloration basée sur densité (échelle log)

### Détection de cycle périodique
- Méthode de Floyd (tortoise and hare)
- Méthode de Newton pour localiser cycles exacts
- Approximation après un nombre suffisant d'itérations

### Paramètres par défaut recommandés
- **Domaine** : `[-2.5, 1.5] × [-1.5, 1.5]` (domaine Mandelbrot)
- **iterationMax** : 10000-100000 (beaucoup plus élevé que escape-time)
- **bailout** : N/A (pas d'échappement)
- **seed** : `(0, 0)`
- **zoomfactor** : 8
- **Échantillonnage** : Aléatoire ou grille dans le Mandelbrot set
- **Détection cycle** : Après 50-100 itérations de warmup

### Références
- **Source originale** : Complémentaire au Buddhabrot (Melinda Green)
- **Code de référence** :
  - Mathr.co.uk : Blog sur Ultimate Anti-Buddhabrot
  - Paul Bourke : Articles sur Buddhabrot
- **Documentation** :
  - Mathr.co.uk : "Ultimate Anti-Buddhabrot"
  - Articles sur mesure invariante et cycles attractifs

### Notes techniques
- **GMP** : Oui (pour précision dans détection de cycles)
- **SIMD** : Partiel (algorithmes spécialisés)
- **Particularités** :
  - Algorithme complètement différent de escape-time
  - Nécessite détection de cycles périodiques
  - Densité accumulée sur cycles attractifs
  - Complémentaire au Buddhabrot (type 16)
  - Nécessite modifications architecturales importantes (comme Buddhabrot)

### Images de référence
- Formes "négatives" du Buddhabrot
- Accumulation de densité sur cycles périodiques
- Structures différentes de l'escape-time standard

---

## 9. Orbit Trap Fractals (Généralisés)

### Concept
Système généralisé d'orbit traps permettant différents types de pièges géométriques :
- **Cross trap** : Axes Re=0 et Im=0 (Pickover stalks)
- **Circle trap** : Cercle centré à l'origine
- **Line trap** : Ligne droite
- **Polygon trap** : Polygone régulier
- **Multiple traps** : Combinaison de plusieurs traps

### Formule de base
```
z₀ = 0
z_{n+1} = z² + c  (ou autre formule)
```

Pendant l'itération, calculer distance à chaque trap :
```
distance_trap_i = distance(z, trap_i)
trap_min = min(trap_min, min_i(distance_trap_i))
```

### Coloration
Utiliser `trap_min` ou combinaison de distances pour coloration.
Peut combiner avec escape-time pour effets mixtes.

### Paramètres par défaut recommandés
- **Domaine** : `[-2.0, 1.0] × [-1.5, 1.5]`
- **iterationMax** : 100-1000
- **bailout** : 4.0
- **Type de trap** : Configurable (cross, circle, line, polygon)
- **Paramètres trap** : 
  - Circle : rayon
  - Line : équation de droite
  - Polygon : nombre de côtés, rayon
- **zoomfactor** : 4-8

### Références
- **Source originale** : Généralisation de Pickover stalks
- **Code de référence** :
  - UltraFractal : Système d'orbit traps
  - FractalForums : Discussions sur traps multiples
- **Documentation** :
  - UltraFractal : Documentation orbit traps
  - FractalForums : Formules et exemples

### Notes techniques
- **GMP** : Oui (si nécessaire)
- **SIMD** : Oui (calculs de distance)
- **Particularités** :
  - Système flexible et extensible
  - Permet grande variété de formes
  - Compatible avec escape-time mais nécessite système de traps
  - Coloration basée sur distances plutôt que iterations
  - Peut être combiné avec escape-time pour effets mixtes

### Images de référence
- Formes géométriques imbriquées
- Patterns variés selon type de trap
- Structures complexes avec détails fins

---

## Résumé des Priorités d'Implémentation

### Priorité 1 (Facile, implémentation immédiate)
1. **Perpendicular Burning Ship** - ⭐ Faible complexité
2. **Celtic Fractal** - ⭐ Faible complexité
3. **Alpha Mandelbrot** - ⭐ Faible complexité

### Priorité 2 (Moyenne complexité, intéressantes)
4. **Pickover Stalks** - ⭐⭐ Moyenne (système orbit trap)
5. **Nova Fractal** - ⭐⭐ Moyenne (dérivée polynomiale)
6. **Multibrot non-entier** - ⭐⭐ Moyenne (branch cuts)

### Priorité 3 (Évaluation approfondie)
7. **Lambda Fractal** - ⭐ Faible mais nécessite validation
8. **Anti-Buddhabrot** - ⭐⭐ Moyenne (algorithme spécialisé)
9. **Orbit Trap généralisés** - ⭐⭐ Moyenne (système flexible)

---

## Notes Générales

### Compatibilité Architecture
Toutes les fractales listées sont compatibles avec l'architecture escape-time de fractall, sauf :
- **Anti-Buddhabrot** : Nécessite algorithme spécialisé (comme Buddhabrot type 16)
- **Pickover Stalks / Orbit Traps** : Nécessitent système de traps mais peuvent utiliser escape-time comme base

### Support GMP
Toutes les fractales peuvent bénéficier du support GMP pour zooms profonds, particulièrement :
- Multibrot (précision importante pour branch cuts)
- Nova (précision dans calculs de dérivée)
- Anti-Buddhabrot (détection précise de cycles)

### Optimisations SIMD
La plupart des fractales peuvent bénéficier d'optimisations SIMD :
- Perpendicular Burning Ship : Oui
- Celtic : Oui
- Alpha Mandelbrot : Oui
- Multibrot : Partiel (exponentiation complexe)
- Nova : Partiel (divisions complexes)

---

## Prochaines Étapes

1. **Phase 2** : Prototypage Python pour valider les formules
2. **Phase 3** : Implémentation C dans fractall
3. **Tests** : Validation avec toutes les palettes et différents paramètres
