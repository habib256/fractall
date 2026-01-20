# Bugs identifiés dans fractall

**Date** : $(date)  
**Statut** : Tous les bugs critiques et mineurs ont été corrigés ✅

## Bugs critiques (corrigés)

### 1. **Bug main.c ligne 375** : Barre d'état non mise à jour après changement de colorRepeat
**Fichier** : `src/main.c` ligne 375  
**Problème** : Après avoir changé le nombre de répétitions du gradient avec la touche R, la barre d'état n'est pas mise à jour pour afficher les nouvelles informations.  
**Impact** : L'utilisateur ne voit pas la mise à jour de la barre d'état après avoir appuyé sur R.  
**Solution** : Ajouter l'appel à `SDLGUI_StateBar_Update` après le changement de colorRepeat, similaire à ce qui est fait pour la touche C.

### 2. **Bug main.c lignes 749, 996, 1086** : Accès potentiel à zmatrix_gmp[0] sans vérification
**Fichier** : `src/main.c` lignes 749, 996, 1086  
**Problème** : Le code accède à `f->zmatrix_gmp[0].x` pour lire la précision sans vérifier que le tableau n'est pas vide (`f->xpixel * f->ypixel > 0` est vérifié, mais pas que l'élément [0] est initialisé).  
**Impact** : Risque de segfault si zmatrix_gmp est alloué mais non initialisé.  
**Solution** : Vérifier que zmatrix_gmp[0] est initialisé avant d'accéder à sa précision, ou utiliser une variable de suivi de précision.

### 3. **Bug precision_detector.c ligne 32** : Pas de vérification division par zéro
**Fichier** : `src/precision_detector.c` ligne 32  
**Problème** : Le code divise par `f->xpixel` et `f->ypixel` sans vérifier qu'ils ne sont pas 0.  
**Impact** : Division par zéro si xpixel ou ypixel est 0, causant un crash.  
**Solution** : Ajouter une vérification que xpixel et ypixel sont > 0 avant de diviser.

## Bugs mineurs (corrigés)

### 4. **Bug main.c ligne 114** : Message d'aide incorrect
**Fichier** : `src/main.c` ligne 114  
**Problème** : Le message d'aide indique que la hauteur GUI par défaut est 66, mais le code utilise 51.  
**Impact** : Information incorrecte pour l'utilisateur.  
**Solution** : Corriger le message pour indiquer 51, ou changer la valeur par défaut à 66.

### 5. **Bug main.c ligne 135** : stateH non recalculé quand menuH change
**Fichier** : `src/main.c` ligne 135  
**Problème** : Quand l'utilisateur change la hauteur du menu avec `-g`, `stateH` n'est pas recalculé.  
**Impact** : La barre d'état peut être mal positionnée si menuH change après l'initialisation.  
**Solution** : Recalculer stateH après avoir modifié menuH.

### 6. **Bug SDLGUI.c ligne 40** : Division entière peut donner 0 pour barw
**Fichier** : `src/SDLGUI.c` ligne 40  
**Problème** : Le calcul `g.barw = g.w*(g.w-4)/((g.buttonSize+2)*buttonNumber)` utilise une division entière qui peut donner 0 si le dénominateur est grand.  
**Impact** : La barre de défilement peut avoir une largeur de 0, rendant le GUI inutilisable.  
**Solution** : Utiliser une division en virgule flottante et arrondir, ou ajouter une valeur minimale.

### 7. **Bug SDLGUI.c ligne 35** : xplaceMax peut être négatif
**Fichier** : `src/SDLGUI.c` ligne 35  
**Problème** : Le calcul `g.xplaceMax = (2 + (g.buttonSize + 2) * g.buttonNumber) - g.w` peut donner une valeur négative si `g.w` est grand et `buttonNumber` petit.  
**Impact** : Comportement imprévisible dans le calcul de la position de la barre de défilement.  
**Solution** : Ajouter une vérification pour s'assurer que xplaceMax >= 0.

## Bug critique découvert lors du test (corrigé)

### 9. **Bug SDLGUI.c ligne 131** : Division par zéro dans SDLGUI_MenuBar_Draw
**Fichier** : `src/SDLGUI.c` ligne 131  
**Problème** : Division par `g->xplaceMax` qui peut être 0 si tous les boutons rentrent dans l'écran, causant une exception en point flottant au démarrage.  
**Impact** : Crash immédiat au démarrage avec "Exception en point flottant (core dumped)".  
**Solution** : Ajout d'une vérification pour éviter la division si `xplaceMax == 0`. Protection supplémentaire pour `buttonSize` négatif et dans `SDLGUI_ButtonNumberRead`.

## Notes

- Les coordonnées GMP sont toujours initialisées dans `Fractal_ChangeType`, donc la libération dans `Fractal_Destroy` est correcte.
- Les vérifications de `f->xpixel * f->ypixel > 0` sont présentes dans la plupart des endroits critiques, mais pas partout.
- **Tous les bugs ont été corrigés et le code compile correctement.**
