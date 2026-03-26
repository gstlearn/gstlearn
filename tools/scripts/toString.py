#!/usr/bin/env python3
import re
import sys
import os

# The aim of this script is to indentify all the classes which inherits
# (recursively) from AStringable in order to add the extension
#
# %extend {class_name} {{
#  std::string __repr__() {{
#    return $self->toString();
#  }}
# 
# in the toString.i file such that toString is automatically called in the target
# languages when the name of the object is given.

def extract_class_and_bases(all_lines):
    """
    Analyzes all lines of code and extracts the class name and the list of base classes inherited,
    while removing the word 'public'. It handles all cases: mono or multiline class declaration.
    The class must be exported using GSTLEARN_EXPORT.

    :param: all_lines: List of lines containing the class header
    :return: Tuple (class name, list of inherited base classes)
    """
    # Regex to extract the class name and all base classes (handle multiline inheritance list)
    pattern = r"class\s+GSTLEARN_EXPORT\s+(\w+)\s*(?::\s*([^{]+))?"
    # Catch all lines (using DOTALL)
    match = re.search(pattern, all_lines, re.DOTALL)

    if match:
        # Class name
        class_name = match.group(1)

        # Inherited base classes, separated by commas and or carriage returns (if any)
        base_classes_raw = match.group(2)
        if base_classes_raw:
            base_classes_raw = base_classes_raw.split(",")
            # Remove the word 'public' and 'virtual' from each base class
            base_classes = [
                base.replace("public", "").replace("virtual", "").strip()
                for base in base_classes_raw
            ]
        else:
            base_classes = []

        return class_name, base_classes
    else:
        return None, []


def get_all_descendants(target_class, hierarchy):
    # Init list of descendants found
    descendants = set()

    # Parse the dictionary and expand the parents list
    # until we find the target class or exhaust the hierarchy
    for child, parents in hierarchy.items():
        if len(parents) == 0:
            continue
        if target_class in parents:
            descendants.add(child)
            continue
        while (target_class not in parents) and len(parents) > 0:
            new_parents = set()
            for parent in parents:
                if parent in hierarchy:
                    new_parents.update(hierarchy[parent])
            if (parents == new_parents) or len(new_parents) == 0:
                break  # No new parents found, stop the search
            parents = new_parents
            if target_class in parents:
                descendants.add(child)
                break  # Target class found, stop the search for this child

    return list(descendants)


def find_classes_inheriting_from_AStringable(root_folder):
    """
    Finds classes that inherit directly from AStringable in all .hpp files within a directory.

    :param root_folder: Folder containing .hpp files (e.g., "include").
    :return: List of classes inheriting directly from AStringable.
    """
    direct_inheritors = []

    # Dictionary to store class hierarchies
    class_hierarchy = {}

    # Traverse through the directory tree
    for root, _, files in os.walk(root_folder):
        for file in files:
            if file.endswith(".hpp"):
                file_path = os.path.join(root, file)

                try:
                    with open(file_path, "r", encoding="utf-8") as f:
                        lines = f.readlines()
                        # Extract class and base classes
                        class_name, base_classes = extract_class_and_bases(
                            "".join(lines)
                        )
                        if class_name:
                            class_hierarchy[class_name] = base_classes
                            # Check if the class inherits directly from AStringable
                            if "AStringable" in base_classes:
                                direct_inheritors.append(class_name)
                except (UnicodeDecodeError, FileNotFoundError):
                    continue

    all_inheritors = get_all_descendants("AStringable", class_hierarchy)

    # Return the classes sorted alphabetically
    return all_inheritors


# Display all classes inheriting from AStringable (directly or indirectly) in alphabetical order


def generate_swig_extend_code(class_name):
    """
    Write SWIG code %extend for a given class.

    :param class_name: Name of the class.
    :return: Code SWIG %extend.
    """
    return f"""
%extend gstlrn::{class_name} {{
  std::string __repr__() {{
    return $self->toString();
  }}
}}
"""


if __name__ == "__main__":
    output_txt_file = sys.argv[1]
    include_path = os.path.join("..", "..", "include")
    Astringable_classes = sorted(find_classes_inheriting_from_AStringable(include_path))
    with open(output_txt_file, "w", encoding="utf-8") as file:
        for class_name in Astringable_classes:
            file.write(generate_swig_extend_code(class_name))
