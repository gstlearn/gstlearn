import os
import re
def extract_macro_calls(header_file, macro_name):
    # Utilisation de os.path pour manipuler les chemins de manière portable
    classname = os.path.basename(header_file).split(".")[0]
    macro_calls = []
    macro_pattern = rf"^(?!#define\s+{macro_name})\s*{macro_name}\((.*?)\)"
    with open(header_file, 'r') as file:
        for line in file:
            match = re.findall(macro_pattern, line, re.MULTILINE)
            if match:
                contents = match[0].split(",")
                if len(contents) == 2:
                    body, func_name = contents
                    arg = ""
                else:
                    if len(contents) == 3:
                        body, func_name, arg = contents
                    else:
                        print("error.")  
                macro_calls += [(classname, func_name.strip(), body.strip(), arg)]
    return macro_calls
    
def find_function_return_type_in_file(method_name, header_path):
    """
    Retourne le type de retour nettoyé (sans pointeur, référence, const, etc.) 
    pour une méthode dans un fichier header donné.
    """
    if not os.path.isfile(header_path):
        print(f"❌ Fichier non trouvé : {header_path}")
        return None

    with open(header_path, 'r', encoding='utf-8') as f:
        lines = f.readlines()

    buffer = ""
    collecting = False

    for line in lines:
        stripped = line.strip()
        if not collecting and re.search(rf'\b{re.escape(method_name)}\s*\(', stripped):
            collecting = True
        if collecting:
            buffer += " " + stripped
            if "{" in stripped or ";" in stripped:
                break

    if not buffer:
        print(f"❌ Méthode {method_name} introuvable dans {header_path}")
        return None

    # Nettoyage
    buffer = re.sub(r'//.*', '', buffer)
    buffer = re.sub(r'/\*.*?\*/', '', buffer)

    # Extraction du type
    match = re.match(r'(?:virtual\s+)?([\w\s:<>,*&]+?)\s+\**' + re.escape(method_name) + r'\s*\(', buffer.strip())
    if match:
        return_type = match.group(1).strip()
        # Nettoyage : suppression de const, &, *, etc.
        return_type = return_type.replace('const', '')
        return_type = return_type.replace('&', '')
        return_type = return_type.replace('*', '')
        return_type = return_type.strip()
        return return_type

    print(f"❓ Signature non reconnue pour {method_name} dans {header_path}")
    return None



def find_include_folder_in_file(classname, file_path):
    """
    Dans un fichier donné, retrouve la ligne `#include "Folder/Classname.hpp"`
    et retourne `Folder`.
    """
    if not os.path.isfile(file_path):
        print(f"❌ Fichier non trouvé : {file_path}")
        return None

    pattern = re.compile(rf'#include\s+"([^"]*\/{re.escape(classname)}\.hpp)"')

    with open(file_path, "r", encoding="utf-8") as f:
        for line in f:
            match = pattern.search(line)
            if match:
                include_path = match.group(1)  # Ex. "Covariances/ACov.hpp"
                return os.path.dirname(include_path)  # Ex. "Covariances"

    print(f"❌ Ligne d'#include pour {classname}.hpp non trouvée dans {file_path}")
    return None
   

def find_method_signature_in_header(method_name, header_path):
    if not os.path.isfile(header_path):
        return None

    with open(header_path, "r", encoding="utf-8") as f:
        lines = f.readlines()

    signature_lines = []
    collecting = False
    brace_depth = 0

    for line in lines:
        stripped = line.strip()

        # Ignore les commentaires ou lignes vides
        if not stripped or stripped.startswith("//") or stripped.startswith("/*"):
            continue

        # Démarrage de la collecte si le nom de méthode apparaît
        if not collecting and re.search(rf'\b{re.escape(method_name)}\s*\(', stripped):
            collecting = True

        if collecting:
            signature_lines.append(stripped)

            # Compte les accolades pour éviter de couper au mauvais moment
            brace_depth += stripped.count('{') - stripped.count('}')

            # Fin de la signature : si ; ou { sans autre {
            if (brace_depth == 0 and ';' in stripped) or (brace_depth <= 1 and '{' in stripped and '}' in stripped):
                break

    if not signature_lines:
        return None

    # Fusionner, supprimer les corps de fonctions
    full_signature = ' '.join(signature_lines)
    full_signature = re.sub(r'\{[^{}]*\}', '', full_signature)  # Retire corps simple
    full_signature = re.sub(r'\s+', ' ', full_signature).strip()

    # Force ; si terminé par ') const' ou 'override', etc.
    if not full_signature.endswith(';'):
        full_signature += ';'

    return full_signature
    
    
def extract_base_class_name(lines):
    """
    Extrait le nom de la classe de base d'après une ligne comme :
    class CovList : public CovBase {
    """
    pattern = re.compile(r'class\s+\w+\s*:\s*public\s+(\w+)')
    for line in lines:
        match = pattern.search(line)
        if match:
            return match.group(1)
    return None
    
def extract_base_classes(header_lines):
    """
    Extrait la liste des classes bases à partir des lignes du header.
    Supporte plusieurs bases, format classique : class Derived : public Base1, public Base2 { ...
    """
    bases = []
    class_decl_pattern = re.compile(r'class\s+\w+\s*:\s*([^\\{]+)\{?')
    for line in header_lines:
        match = class_decl_pattern.search(line)
        if match:
            bases_list = match.group(1)
            # extraire les noms de bases (supprimer "public ", "private ", etc.)
            bases = [b.strip().split()[-1] for b in bases_list.split(',')]
            break
    return bases

def find_header_file_recursive(class_name, root_dir):
    """
    Recherche récursive d'un fichier <class_name>.hpp dans root_dir et ses sous-dossiers.
    Renvoie le chemin complet du fichier trouvé ou None.
    """
    target = f"{class_name}.hpp"
    for dirpath, _, files in os.walk(root_dir):
        if target in files:
            return os.path.join(dirpath, target)
    return None
    
def extract_class_and_bases(line):
    pattern = r'class\s+\w+\s+(\w+)\s*[:\s]*([^{\n]+)?\s*(\{)?'
    match = re.match(pattern, line)
    if match:
        class_name = match.group(1)
        base_classes_raw = match.group(2)
        if base_classes_raw:
            base_classes_raw = base_classes_raw.split(',')
            base_classes = [b.replace('public', '').strip() for b in base_classes_raw]
        else:
            base_classes = []
        return class_name, base_classes
    else:
        return None, []

def find_header_for_class(class_name, root_folder):
    expected_file = f"{class_name}.hpp"
    for dirpath, _, filenames in os.walk(root_folder):
        if expected_file in filenames:
            return os.path.join(dirpath, expected_file)
    return None

def find_method_signature_in_header(method_name, header_path):
    if not os.path.isfile(header_path):
        return None

    with open(header_path, "r", encoding="utf-8") as f:
        lines = f.readlines()

    signature_lines = []
    collecting = False
    brace_depth = 0

    for line in lines:
        stripped = line.strip()

        if not stripped or stripped.startswith("//") or stripped.startswith("/*"):
            continue

        if not collecting and re.search(rf'\b{re.escape(method_name)}\s*\(', stripped):
            collecting = True

        if collecting:
            signature_lines.append(stripped)
            brace_depth += stripped.count('{') - stripped.count('}')
            if (brace_depth == 0 and ';' in stripped) or (brace_depth <= 1 and '{' in stripped and '}' in stripped):
                break

    if not signature_lines:
        return None

    full_signature = ' '.join(signature_lines)
    full_signature = re.sub(r'\{[^{}]*\}', '', full_signature)
    full_signature = re.sub(r'\s+', ' ', full_signature).strip()

    if not full_signature.endswith(';'):
        full_signature += ';'

    return full_signature

def find_method_in_class_or_bases(class_name, method_name, root_folder, visited=None, depth=0):
    if visited is None:
        visited = set()
    if class_name in visited:
        return None
    visited.add(class_name)

    header_path = find_header_for_class(class_name, root_folder)
    if not header_path:
        print(f"[INFO] Header not found for class {class_name}")
        return None

    # 1. Cherche la méthode dans le header direct
    signature = find_method_signature_in_header(method_name, header_path)
    if signature:
        if depth > 0:
            print(f"[INFO] Méthode '{method_name}' trouvée dans la base '{class_name}'")
        return signature

    # 2. Si non trouvée, cherche dans les bases
    try:
        with open(header_path, 'r', encoding='utf-8') as f:
            lines = f.readlines()
    except Exception:
        return None

    for line in lines:
        _, base_classes = extract_class_and_bases(line)
        for base_class in base_classes:
            signature = find_method_in_class_or_bases(base_class, method_name, root_folder, visited, depth + 1)
            if signature:
                return signature

    return None   

