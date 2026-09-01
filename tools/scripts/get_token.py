#!/usr/bin/env python3
import urllib.parse
import webbrowser
import os


def main():
    # Configuration des paramètres du jeton
    note = "Script-NonReg-Times"
    scopes = ["repo"]  # Coche automatiquement la case repo

    params = {"description": note, "scopes": ",".join(scopes)}

    # URL de création de jeton pré-remplie
    url = f"https://github.com/settings/tokens/new?{urllib.parse.urlencode(params)}"

    print("Ouverture de la page GitHub de création de jeton...")
    webbrowser.open(url)

    print("\n--- INSTRUCTIONS ---")
    print("1. Valide la création au bas de la page web.")
    print("2. Copie le jeton généré.")

    token = input("\nColle ton jeton ici : ").strip()

    if token:
        print("\nPour l'utiliser dans ta session csh actuelle, exécute :")
        print(f"setenv GITHUB_TOKEN {token}")

        # Optionnel : enregistre la commande dans un fichier à sourcer
        with open("set_env.csh", "w") as f:
            f.write(f"setenv GITHUB_TOKEN {token}\n")
        print("\nUn fichier 'set_env.csh' a été créé. Tu peux faire :")
        print("source set_env.csh")


if __name__ == "__main__":
    main()
