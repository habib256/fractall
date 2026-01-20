# Plan de Diagnostic et Correction : Nova Fractal

## Problèmes identifiés potentiels

### 1. Condition de bailout incorrecte

**Problème actuel** :
- Utilise `Magz(z) < f.bailout` avec `bailout = 100.0`
- C'est une condition d'échappement (comme Mandelbrot)
- Nova est une fractale de **convergence**, pas d'échappement

**Impact** : La fractale ne détecte pas la convergence vers les racines, ce qui peut donner des résultats incorrects.

**Solution** : Implémenter une détection de convergence basée sur :
- Distance entre itérations successives : `|z_{n+1} - z_n| < epsilon`
- Distance normalisée : `|z_{n+1} - z_n| / max(1, |z_n|) < epsilon`
- Distance aux racines connues : `|z_n - root_i| < tolerance`

### 2. Valeur initiale z0

**Problème actuel** :
- `z = zPixel` (le point dans le plan complexe)
- Selon certaines sources, pour Nova Mandelbrot, `z0` devrait être fixe (ex: `z0 = 1` ou point critique)

**Impact** : La valeur initiale peut affecter significativement le résultat.

**Solutions à tester** :
- Option A : `z0 = 1` (racine de `z^3 - 1 = 0`)
- Option B : `z0 = zPixel` (comme actuellement)
- Option C : `z0 = point critique` (racine de `p'(z) = 0`)

### 3. Coloration inadaptée

**Problème actuel** :
- Coloration basée sur `iteration` (nombre d'itérations jusqu'à échappement)
- Pour Nova, on devrait colorer selon :
  - La racine vers laquelle converge le point
  - La vitesse de convergence
  - La distance à la racine

**Impact** : Les couleurs peuvent ne pas refléter correctement la structure de la fractale.

### 4. Comparaison avec Newton

**Observation** :
- Newton utilise aussi `Magz(z) < f.bailout` avec `bailout = 4`
- Newton fonctionne-t-il correctement ? Si oui, pourquoi Nova ne fonctionne pas ?
- Si Newton ne fonctionne pas non plus, le problème est systémique

## Plan de diagnostic

### Étape 1 : Vérifier la formule mathématique

**Actions** :
1. Vérifier que la formule implémentée correspond exactement à :
   ```
   z_{n+1} = z_n - a·((z_n^p - 1)/(p·z_n^(p-1))) + c
   ```
2. Vérifier l'ordre des opérations :
   - Calcul de `z^p` et `z^(p-1)`
   - Calcul de `(z^p - 1) / (p·z^(p-1))`
   - Multiplication par `a`
   - Soustraction de `z`
   - Addition de `c`

**Fichier** : `src/EscapeTime.c` lignes 2529-2569

**Vérifications** :
- [ ] La formule est mathématiquement correcte
- [ ] L'ordre des opérations est correct
- [ ] Les calculs intermédiaires sont corrects

### Étape 2 : Tester différentes valeurs initiales

**Actions** :
1. Créer des versions de test avec différentes valeurs initiales :
   - Version A : `z = MakeComplex(1.0, 0.0)` (racine de z³-1=0)
   - Version B : `z = zPixel` (actuel)
   - Version C : `z = MakeComplex(0.0, 0.0)` puis skip première itération si nécessaire
2. Comparer les résultats visuels
3. Identifier quelle valeur donne le résultat attendu

**Tests à effectuer** :
- [ ] Tester avec `z0 = 1`
- [ ] Tester avec `z0 = zPixel`
- [ ] Comparer visuellement les résultats
- [ ] Identifier la meilleure option

### Étape 3 : Analyser le comportement de convergence

**Actions** :
1. Ajouter des logs de débogage pour tracer :
   - Valeur de `z` à chaque itération
   - Distance `|z_{n+1} - z_n|`
   - Distance aux racines connues
2. Identifier si la convergence se produit
3. Vérifier si le bailout est atteint trop tôt ou trop tard

**Points à vérifier** :
- [ ] Les points convergent-ils vers les racines ?
- [ ] La distance entre itérations diminue-t-elle ?
- [ ] Le bailout est-il approprié ?

### Étape 4 : Comparer avec implémentations de référence

**Actions** :
1. Rechercher des implémentations de référence (UltraFractal, FractalShades, etc.)
2. Comparer :
   - Formule exacte
   - Valeur initiale
   - Condition de bailout
   - Paramètres par défaut
3. Identifier les différences

**Sources de référence** :
- UltraFractal documentation
- FractalShades code
- HPDZ.net examples
- Mandelbrowser tutorials

### Étape 5 : Tester différents paramètres

**Actions** :
1. Tester avec `a = 1.0` (actuel)
2. Tester avec `a = 0.5` (convergence plus lente)
3. Tester avec `a = 2.0` (convergence plus rapide, peut créer chaos)
4. Tester avec `p = 2` au lieu de `p = 3`
5. Tester avec différents `bailout` (4.0, 10.0, 100.0)

**Paramètres à tester** :
- [ ] `a = 0.5, 1.0, 2.0`
- [ ] `p = 2, 3, 4`
- [ ] `bailout = 4, 10, 20, 100`
- [ ] `iterationMax = 100, 500, 1000`

## Corrections proposées

### Correction 1 : Implémenter détection de convergence

**Fichier** : `src/EscapeTime.c`

**Modification** : Remplacer la condition de bailout par une détection de convergence

```c
// AVANT :
} while ((i < f.iterationMax) && ( Magz (z) < f.bailout));

// APRÈS :
complex z_prev;
double conv_epsilon = 1e-7;  // Seuil de convergence
z_prev = z;
do {
    // ... calcul de z ...
    
    // Détection de convergence
    double diff_sq = Magz2(Subz(z, z_prev));
    double z_sq = Magz2(z);
    double denom = (z_sq < 1.0) ? 1.0 : z_sq;
    
    if (diff_sq / denom < conv_epsilon * conv_epsilon) {
        // Convergence détectée
        break;
    }
    
    z_prev = z;
    
    // Échappement si |z| devient trop grand
    if (Magz(z) > f.bailout) {
        break;
    }
} while (i < f.iterationMax);
```

### Correction 2 : Ajuster valeur initiale

**Option recommandée** : Utiliser `z0 = 1` pour Nova Mandelbrot

```c
// Pour Nova Mandelbrot, utiliser z0 = 1 (racine de z^3 - 1 = 0)
z = MakeComplex(1.0, 0.0);
```

**Alternative** : Utiliser `z0 = zPixel` mais avec gestion spéciale du premier appel.

### Correction 3 : Ajuster paramètres par défaut

**Modifications** :
- `bailout` : Réduire à `10.0` ou `20.0` (au lieu de 100.0)
- `iterationMax` : Augmenter à `500` ou `1000` pour capturer plus de détails
- Ajouter un paramètre `conv_epsilon` pour la tolérance de convergence

### Correction 4 : Améliorer la coloration (optionnel)

**Modification** : Stocker la racine vers laquelle converge le point dans `result.iteration` ou utiliser `result.z` pour déterminer la racine.

## Tests de validation

1. **Test visuel** :
   - Comparer avec des images de référence de Nova fractals
   - Vérifier la présence de structures en spirales
   - Vérifier la présence de "bras" ou "antennes"

2. **Test de convergence** :
   - Vérifier que les points convergent vers les racines attendues
   - Vérifier que la détection de convergence fonctionne

3. **Test de performance** :
   - Vérifier que le temps de rendu reste acceptable
   - Comparer avec d'autres fractales similaires (Newton)

## Ordre d'implémentation recommandé

1. **Priorité 1** : Corriger la valeur initiale (`z0 = 1`)
2. **Priorité 2** : Implémenter détection de convergence
3. **Priorité 3** : Ajuster les paramètres par défaut
4. **Priorité 4** : Améliorer la coloration (si nécessaire)

## Notes techniques

- La fonction `Magz2()` existe déjà dans `complexmath.h` pour calculer `|z|²` efficacement
- La détection de convergence nécessite de stocker `z_prev` à chaque itération
- Pour `p=3`, les racines de `z³ - 1 = 0` sont : `1`, `(-1 ± i√3)/2`
- Le système de colorisation actuel peut être adapté pour utiliser la racine de convergence

## Hypothèses à vérifier

1. **Hypothèse 1** : La valeur initiale `z0 = zPixel` n'est pas appropriée
   - **Test** : Changer pour `z0 = 1`
   - **Résultat attendu** : Fractale plus proche des références

2. **Hypothèse 2** : Le bailout de 100.0 est trop élevé
   - **Test** : Réduire à 10.0 ou 20.0
   - **Résultat attendu** : Meilleure détection des zones intéressantes

3. **Hypothèse 3** : La détection de convergence est nécessaire
   - **Test** : Implémenter détection de convergence
   - **Résultat attendu** : Fractale plus précise et détaillée

4. **Hypothèse 4** : Le paramètre `a` ou `p` doit être ajusté
   - **Test** : Tester différentes valeurs
   - **Résultat attendu** : Structures plus visibles

## Références

- Formule standard : `z_{n+1} = z_n - a·((z_n^p - 1)/(p·z_n^(p-1))) + c`
- Pour `p=3` : `z_{n+1} = z_n - a·((z_n³ - 1)/(3·z_n²)) + c`
- Racines de `z³ - 1 = 0` : `1`, `(-1 ± i√3)/2`
- Valeur initiale recommandée : `z0 = 1` pour Nova Mandelbrot
