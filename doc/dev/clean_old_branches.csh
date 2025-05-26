# This script is meant to delete the VSCODE branches that:
# - you may still have on your computer
# - that have already been merged into 'dev'
# - that do not have any new material since the last merge
# These branches can be safely deleted

set DIR = $GSTLEARN_DIR/gstlearn
cd $DIR

# Récupère la liste des branches (avec * éventuellement)
set list = `git branch --format="%(refname:short)"`

# Nettoie les noms de branches (supprime '* ' éventuel)
set clean_list = ()
foreach b ($list)
    set b = `echo $b | sed 's/^\* //'`
    set clean_list = ($clean_list $b)
end

# Vérifie que la branche 'dev' existe
set dev_exists = 0
foreach b ($clean_list)
    if ("$b" == "dev") then
        set dev_exists = 1
        break
    endif
end

if (! $dev_exists) then
    echo "Erreur : la branche 'dev' n'existe pas dans ce dépôt local."
    exit 1
endif

# Boucle sur toutes les branches
foreach branch ($clean_list)

    # On ne propose pas la branche 'dev' elle-même
    if ("$branch" == "dev") continue

    # On ne propose pas la branche courante (protège contre autodestruction)
    set current = `git branch --show-current`
    if ("$branch" == "$current") continue

    # Vérifie si 'dev' contient cette branche
    set contains_list = `git branch --contains $branch`
    set dev_found = 0
    foreach b ($contains_list)
        set b = `echo $b | sed 's/^\* //'`
        if ("$b" == "dev") then
            set dev_found = 1
            break
        endif
    end

    if ($dev_found) then
        echo " "
        echo "La branche '$branch' peut être détruire."
        echo -n "Souhaitez-vous la supprimer ? (o/n[def]) : "
        set confirmation = $<
        if ("$confirmation" == "o" || "$confirmation" == "O") then
            git branch -d "$branch"
        else
            echo "Branche '$branch' conservée."
        endif
    else
        echo "Branche '$branch' n'a pas été mergée à 'dev'."
    endif
end
