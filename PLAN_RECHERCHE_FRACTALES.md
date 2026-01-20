# Plan de Recherche : Nouveaux Types de Fractales pour fractall

## 📋 Vue d'ensemble

Ce document présente un plan structuré pour identifier, rechercher et implémenter de nouveaux types de fractales visuellement intéressantes pour le projet fractall.

**État actuel** : 17 types de fractales (1-2 vectorielles, 3-17 escape-time)

**Objectif** : Identifier 5-10 nouveaux types de fractales esthétiquement remarquables et techniquement réalisables avec l'architecture existante.

---

## 🎯 Critères de sélection

### Critères visuels
- ✅ Formes distinctes et reconnaissables
- ✅ Potentiel de zoom profond avec détails intéressants
- ✅ Compatibilité avec les 8 palettes existantes
- ✅ Esthétique variée (organique, géométrique, abstraite)

### Critères techniques
- ✅ Compatible avec l'architecture escape-time existante
- ✅ Calculable avec double/GMP (précision arbitraire)
- ✅ Optimisable avec SIMD/OpenMP
- ✅ Formule mathématique claire et documentée
- ✅ Paramètres ajustables (seed, bailout, etc.)

### Critères d'implémentation
- ✅ Complexité modérée (pas de dépendances externes lourdes)
- ✅ Performance acceptable (< 10s pour rendu initial à 800x600)
- ✅ Compatible avec Divergence Detection (DDp1/DDp2)

---

## 🔍 Catégories de recherche

### Catégorie 1 : Variantes Mandelbrot/Julia (Priorité HAUTE)

Ces fractales sont des variations directes des formules classiques, donc faciles à intégrer.

#### 1.1 Perpendicular Burning Ship
- **Formule** : `z(n+1) = (Re(z) - i×|Im(z)|)² + c`
- **Pourquoi** : Variante du Burning Ship (type 13) avec pliage perpendiculaire
- **Complexité** : ⭐ Faible (similaire à Burning Ship)
- **Esthétique** : Formes géométriques nettes, structures symétriques
- **Sources** : GitHub kf2, FractalShades
- **Statut** : 🟢 Recherche complétée - Voir FRACTALES_RECHERCHE.md

#### 1.2 Celtic Fractal
- **Formule** : `z(n+1) = |Re(z²)| + i×|Im(z²)| + c` (variante)
- **Pourquoi** : Formes celtiques distinctes, très populaire
- **Complexité** : ⭐ Faible
- **Esthétique** : Motifs entrelacés, symétrie complexe
- **Sources** : DeviantArt, UltraFractal
- **Statut** : 🟢 Recherche complétée - Voir FRACTALES_RECHERCHE.md

#### 1.3 Multibrot (puissances non-entières)
- **Formule** : `z(n+1) = z(n)^d + c` où d est réel (2.5, 3.7, etc.)
- **Pourquoi** : Morphing entre formes, animations possibles
- **Complexité** : ⭐⭐ Moyenne (gestion branch cuts)
- **Esthétique** : Transitions fluides, symétries variables
- **Sources** : Paul Bourke, Mitch Reid
- **Statut** : 🟢 Recherche complétée - Voir FRACTALES_RECHERCHE.md

#### 1.4 Alpha Mandelbrot (nested/composite)
- **Formule** : `z(n+1) = z² + (z² + c)² + c`
- **Pourquoi** : Structures superposées, détails supplémentaires
- **Complexité** : ⭐ Faible
- **Esthétique** : Multi-couches, auto-similarité renforcée
- **Sources** : Reddit r/fractals
- **Statut** : 🟢 Recherche complétée - Voir FRACTALES_RECHERCHE.md

---

### Catégorie 2 : Fractales basées sur fonctions spéciales (Priorité MOYENNE)

Utilisation de fonctions transcendantes (sin, exp, log) pour créer des formes organiques.

#### 2.1 Pickover Stalks / Biomorphs
- **Formule** : `z(n+1) = sin(z) + exp(z) + c` ou variantes
- **Coloration** : Orbit trap sur axes (min(|Re|, |Im|))
- **Pourquoi** : Formes biologiques/organiques uniques
- **Complexité** : ⭐⭐ Moyenne (orbit trap spécial)
- **Esthétique** : Formes biomorphiques, structures en croix
- **Sources** : Paul Bourke, Wikipedia
- **Statut** : 🟢 Recherche complétée - Voir FRACTALES_RECHERCHE.md

#### 2.2 Nova Fractal (variante Newton)
- **Formule** : `z(n+1) = z - a×p(z)/p'(z) + c` où p est polynôme
- **Pourquoi** : Spirales élégantes, structures en antennes
- **Complexité** : ⭐⭐ Moyenne (dérivée polynomiale)
- **Esthétique** : Bras spiralés, motifs répétitifs fins
- **Sources** : UltraFractal, HPDZ.net
- **Statut** : 🟢 Recherche complétée - Voir FRACTALES_RECHERCHE.md

#### 2.3 Lambda Fractal
- **Formule** : `z(n+1) = λ×z×(1-z)` (suite logistique complexe)
- **Pourquoi** : Complémentaire à Lyapunov (type 17)
- **Complexité** : ⭐ Faible
- **Esthétique** : Formes organiques, connexions avec chaos
- **Sources** : Wikipedia, mathématiques dynamiques
- **Statut** : 🟢 Recherche complétée - Voir FRACTALES_RECHERCHE.md

---

### Catégorie 3 : Fractales 3D/Quaternions (Priorité BASSE)

Ces fractales nécessitent des modifications architecturales plus importantes.

#### 3.1 Mandelbox
- **Formule** : Opérations de pliage/échelle sur quaternions
- **Pourquoi** : Structures 3D fascinantes
- **Complexité** : ⭐⭐⭐ Élevée (rendu 3D requis)
- **Esthétique** : Boîtes récursives, structures 3D
- **Sources** : FractalForums, Tom Lowe
- **Statut** : 🔴 À évaluer (peut nécessiter rendu 3D)

#### 3.2 Quaternion Julia/Mandelbrot
- **Formule** : Extension à 4D (quaternions)
- **Pourquoi** : Exploration dimensionnelle
- **Complexité** : ⭐⭐⭐ Élevée (projection 3D/4D)
- **Esthétique** : Formes 3D complexes
- **Sources** : FractalForums
- **Statut** : 🔴 À évaluer (complexité architecturale)

---

### Catégorie 4 : Fractales algorithmiques spéciales (Priorité VARIABLE)

Algorithmes de rendu différents (comme Buddhabrot).

#### 4.1 Anti-Buddhabrot
- **Algorithme** : Trajectoires des points qui ne s'échappent PAS
- **Pourquoi** : Complémentaire au Buddhabrot (type 16)
- **Complexité** : ⭐⭐ Moyenne (algorithme spécialisé)
- **Esthétique** : Formes négatives du Buddhabrot
- **Sources** : FractalForums, Wikipedia
- **Statut** : 🟢 Recherche complétée - Voir FRACTALES_RECHERCHE.md

#### 4.2 Orbit Trap Fractals (généralisés)
- **Algorithme** : Pièges géométriques multiples (cercles, lignes, polygones)
- **Pourquoi** : Grande variété de formes selon le trap
- **Complexité** : ⭐⭐ Moyenne (système de traps)
- **Esthétique** : Formes géométriques imbriquées
- **Sources** : UltraFractal, FractalForums
- **Statut** : 🟢 Recherche complétée - Voir FRACTALES_RECHERCHE.md

---

## 📊 Plan de recherche par phase

### Phase 1 : Recherche bibliographique (Semaine 1)

#### Tâches
1. **Recherche académique**
   - [ ] Articles scientifiques sur nouvelles fractales (arXiv, Google Scholar)
   - [ ] Thèses sur fractales mathématiques
   - [ ] Livres de référence (Mandelbrot, Peitgen)

2. **Recherche communautaire**
   - [ ] FractalForums.com (base de données de formules)
   - [ ] UltraFractal.com (formules et exemples)
   - [ ] DeviantArt (galeries fractales avec formules)
   - [ ] Reddit r/fractals (discussions récentes)

3. **Recherche code open-source**
   - [ ] GitHub : recherche "fractal formula" "mandelbrot variant"
   - [ ] FractalShades (Python)
   - [ ] Kalles Fraktaler 2 (C++)
   - [ ] Mandelbrot Explorer (JavaScript)

#### Livrables Phase 1
- Liste de 20-30 formules candidates avec références
- Classification par complexité d'implémentation
- Galerie d'images de référence (si disponibles)

---

### Phase 2 : Analyse technique (Semaine 2)

#### Tâches
1. **Évaluation de chaque candidat**
   - [ ] Vérifier compatibilité avec architecture escape-time
   - [ ] Estimer complexité d'implémentation (1-5)
   - [ ] Tester formules mathématiques (prototype Python)
   - [ ] Évaluer performance théorique

2. **Prototypage rapide**
   - [ ] Implémenter 5-10 formules en Python/NumPy
   - [ ] Générer images de test (800x600)
   - [ ] Comparer esthétique et performance
   - [ ] Documenter paramètres intéressants

3. **Sélection finale**
   - [ ] Choisir top 5-10 selon critères
   - [ ] Prioriser par facilité d'implémentation
   - [ ] Vérifier originalité (pas de doublons avec types existants)

#### Livrables Phase 2
- Prototypes Python fonctionnels
- Images de référence pour chaque candidat sélectionné
- Document technique avec formules exactes et paramètres
- Plan d'implémentation par priorité

---

### Phase 3 : Implémentation (Semaines 3-4)

#### Tâches par fractale sélectionnée
1. **Implémentation C**
   - [ ] Ajouter fonction `*_Iteration()` dans EscapeTime.c
   - [ ] Ajouter fonction `*_def()` pour paramètres par défaut
   - [ ] Intégrer dans `FormulaSelector()`
   - [ ] Ajouter nom dans `Fractal_GetTypeName()`

2. **Support GMP** (si nécessaire)
   - [ ] Implémenter version GMP de l'itération
   - [ ] Tester précision arbitraire

3. **Optimisations**
   - [ ] Vérifier compatibilité SIMD (si applicable)
   - [ ] Tester parallélisation OpenMP
   - [ ] Optimiser bailout et iterationMax

4. **Tests et validation**
   - [ ] Rendu à différentes résolutions
   - [ ] Tests de zoom profond
   - [ ] Validation avec toutes les palettes
   - [ ] Comparaison avec références

#### Livrables Phase 3
- Code C intégré et testé
- Documentation dans CLAUDE.md
- Images de démonstration
- Tests de performance

---

## 🎨 Ressources de recherche identifiées

### Sites web de référence
1. **FractalForums.com**
   - Base de données de formules
   - Discussions techniques
   - Galeries d'images

2. **UltraFractal.com**
   - Bibliothèque de formules
   - Documentation détaillée
   - Exemples visuels

3. **Paul Bourke (paulbourke.net)**
   - Articles sur fractales
   - Formules et algorithmes
   - Exemples de code

4. **Wikipedia**
   - Articles sur types de fractales
   - Formules mathématiques
   - Historique et théorie

### Repositories GitHub
1. **smurfix/kf2** (Kalles Fraktaler 2)
   - Implémentations C++ de nombreuses variantes
   - Formules documentées

2. **gbillotey/Fractalshades**
   - Python avec formules variées
   - Documentation technique

3. **simplesummit/fractalexplorer**
   - Explorateur interactif
   - Formules JavaScript

### Articles scientifiques
1. **arXiv.org**
   - Recherche "fractal" + "complex dynamics"
   - Articles récents (2020-2024)

2. **Google Scholar**
   - Recherche "new fractal types"
   - Citations de Mandelbrot, Peitgen

---

## 📝 Template de documentation pour chaque fractale

Pour chaque nouvelle fractale identifiée, documenter :

```markdown
### [Nom de la fractale] (Type X)

**Formule** : `z(n+1) = ...`

**Description** : ...

**Paramètres par défaut** :
- Domaine : [xmin, xmax] × [ymin, ymax]
- iterationMax : ...
- bailout : ...
- seed (si Julia) : ...

**Caractéristiques visuelles** :
- Formes principales : ...
- Zones intéressantes à explorer : ...
- Compatibilité palettes : ...

**Références** :
- Source originale : ...
- Implémentations de référence : ...
- Articles/documentation : ...

**Complexité d'implémentation** : ⭐/⭐⭐⭐⭐⭐

**Notes techniques** :
- Support GMP : Oui/Non
- Optimisations SIMD : Oui/Non
- Particularités : ...
```

---

## 🎯 Priorités recommandées

### Priorité 1 (Implémentation immédiate)
1. **Perpendicular Burning Ship** - Facile, complémentaire au Burning Ship
2. **Celtic Fractal** - Facile, très populaire
3. **Alpha Mandelbrot** - Facile, structures intéressantes

### Priorité 2 (Implémentation à court terme)
4. **Pickover Stalks** - Moyenne, formes uniques
5. **Nova Fractal** - Moyenne, esthétique élégante
6. **Multibrot (non-entier)** - Moyenne, morphing intéressant

### Priorité 3 (Évaluation approfondie)
7. **Lambda Fractal** - Facile mais nécessite validation
8. **Anti-Buddhabrot** - Moyenne, complémentaire au Buddhabrot
9. **Orbit Trap généralisés** - Moyenne, système flexible

### Priorité 4 (Évaluation future)
10. **Fractales 3D/Quaternions** - Complexe, nécessite architecture 3D

---

## ✅ Checklist de validation

Avant d'ajouter une nouvelle fractale au projet :

- [ ] Formule mathématique vérifiée et documentée
- [ ] Prototype Python fonctionnel avec images de référence
- [ ] Implémentation C testée et validée
- [ ] Support GMP (si nécessaire pour zoom profond)
- [ ] Compatible avec toutes les palettes (testé)
- [ ] Performance acceptable (< 10s rendu initial)
- [ ] Documentation ajoutée dans CLAUDE.md
- [ ] Nom ajouté dans Fractal_GetTypeName()
- [ ] Bouton GUI ajouté (si nécessaire)
- [ ] Tests de zoom profond réussis
- [ ] Code review et optimisation

---

## 📚 Références bibliographiques

### Livres
- Mandelbrot, B. B. (1982). *The Fractal Geometry of Nature*
- Peitgen, H.-O., & Richter, P. H. (1986). *The Beauty of Fractals*
- Devaney, R. L. (1992). *A First Course in Chaotic Dynamical Systems*

### Articles clés
- [À compléter lors de la recherche bibliographique]

### Sites web
- FractalForums.com
- UltraFractal.com
- Paul Bourke's Fractals (paulbourke.net/fractals/)
- Wikipedia Fractal articles

---

## 🔄 Mise à jour du plan

Ce plan doit être mis à jour régulièrement :
- Après chaque phase de recherche
- Lors de la découverte de nouvelles ressources
- Après validation/invalidation de candidats
- Lors de l'ajout de nouvelles fractales au projet

**Dernière mise à jour** : [Date]

---

## 📧 Contact et contributions

Pour suggestions ou contributions à ce plan de recherche :
- Ouvrir une issue sur le repository
- Proposer des modifications via pull request
- Partager des ressources intéressantes
