#!/bin/csh

# === À PERSONNALISER ===
set old_branch = "Modif_setLocator"

# === VÉRIFICATION DANS UN RÉPO GIT ===
if (! -d .git) then
    echo "❌ Ce répertoire n'est pas un dépôt Git."
    exit 1
endif

# === MISE À JOUR DE origin/dev ===
echo "📦 Mise à jour de origin/dev..."
git fetch origin dev >& /dev/null

# === VÉRIFICATION DE LA BRANCHE LOCALE ===
git rev-parse --verify $old_branch >& /dev/null
if ($status != 0) then
    echo "❌ Erreur : la branche $old_branch n'existe pas localement."
    exit 1
endif

# === DÉTERMINATION DU MERGE-BASE ===
set base = `git merge-base $old_branch origin/dev | head -n 1`

if ("$base" == "") then
    echo "❌ Impossible de déterminer un merge-base entre $old_branch et origin/dev"
    exit 1
endif

echo "🔍 Merge-base entre $old_branch et origin/dev : $base"

# === GÉNÉRATION DU DIFF ===
set patchfile = "/tmp/diff_from_${old_branch}_to_dev.diff"
git diff $base..$old_branch > $patchfile

# === AFFICHAGE DES RÉSULTATS ===
if ( `wc -l < $patchfile` == 0 ) then
    echo "✅ Aucun changement propre à $old_branch (déjà intégré ou vide)."
    rm -f $patchfile
else
    echo "📄 Modifications propres à $old_branch disponibles dans :"
    echo "   $patchfile"
    echo "💡 Ouvre avec : code $patchfile  ou less $patchfile"
endif
