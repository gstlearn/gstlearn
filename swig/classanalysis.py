#!/usr/bin/env python3
import re
import sys
import os


import re

def find_classes_with_to_string(classname,header_file):
    """
    Analyse un fichier .h pour trouver les classes qui ont une méthode toString.

    :param header_file: Chemin vers le fichier .h.
    :return: Une liste des noms de classes ayant une méthode toString.
    """
    if classname == "String":
        return []
    if classname == "VectorNumT":
        return []
    class_pattern = r"String\s+toString\("
    with open(header_file, 'r', encoding='utf-8') as file:
        for line_number, line in enumerate(file, start=1):
            if re.search(class_pattern, line):
                return [classname]
    return []

def generate_swig_extend_code(class_name):
    """
    Génère le code SWIG %extend pour une classe donnée.

    :param class_name: Nom de la classe.
    :return: Code SWIG %extend.
    """
    return f"""
%extend {class_name} {{
  std::string __repr__() {{
    return $self->toString();
  }}
}}
"""

def add_swig_extensions(classname, header_file, output_file, first = True):
    """
    Analyse un fichier .h, génère le code SWIG %extend pour les classes avec toString,
    et écrit le résultat dans un fichier de sortie.

    :param header_file: Chemin vers le fichier .h.
    :param output_file: Chemin vers le fichier de sortie SWIG.
    """
    classes = find_classes_with_to_string(classname,header_file)
    if not classes:
        return first

    print(f"Classes with toString in {header_file}: {classes}")
    if first:
        mode = 'w'
    else:
        mode = 'a'
    with open(output_file, mode, encoding='utf-8') as file:
        for class_name in classes:
            swig_code = generate_swig_extend_code(class_name)
            file.write(swig_code)
            file.write("\n")  # Ajouter une ligne vide entre les extensions
    return False


def process_cpp_file(classname,header_file, output_file, first):
    pattern = re.compile(r'\w+\s*\*\s*(create\w+)')
    if first:
        mode = 'w'
    else:
        mode = 'a'
    with open(header_file, 'r', encoding='utf-8') as cpp_file, open(output_file, mode, encoding='utf-8') as output:
        for line in cpp_file:
            match = pattern.search(line)
            if match:
                output.write(f"%newobject {classname}::{match.group(1)};\n")
    return False

def extract_included_files(file_path):
    """
    Parcourt un fichier et extrait tous les noms de fichiers inclus avec #include.

    :param file_path: Chemin vers le fichier à analyser.
    :return: Une liste des noms de fichiers inclus.
    """
    included_files = []
    # Expression régulière pour capturer le contenu entre guillemets après #include
    include_pattern = re.compile(r'#include\s+"([^"]+)"')

    with open(file_path, 'r', encoding='utf-8') as file:
        for line in file:
            match = include_pattern.search(line)
            if match:
                # Ajouter le nom du fichier inclus à la liste
                included_files.append(match.group(1))

    return included_files


if __name__ == "__main__":
    filename  = sys.argv[1]
    output_txt_file = sys.argv[2]
    output_txt_file2 = sys.argv[3]
    first = True
    first2 = True
    print("--------------------------------------")
    print(os.getcwd())
    print(output_txt_file)
    print(output_txt_file2)
    files = extract_included_files(filename)
    for file in files:
        fsplit = file.split("/")
        if len(fsplit) != 2:
            continue
        path = os.path.join(fsplit[0], fsplit[1])
        if path[0] in ["include", "Core"]:
            continue
        classname = fsplit[1].split(".")[0]
        header_file = os.path.join("..","..","include", file)  # Utilisation de os.path.join
        first = process_cpp_file(classname,header_file, output_txt_file, first)
        first2 = add_swig_extensions(classname,header_file, output_txt_file2,first2)
