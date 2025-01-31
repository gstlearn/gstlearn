#!/usr/bin/env python3
import re
import sys
import os


def extract_class_and_bases(line):
    """
    Analyzes a line of code and extracts the class name and the list of base classes inherited,
    while removing the word 'public'. It handles both cases: brace on the same line or on the next line.

    :param line: A line containing the class declaration, e.g.:
                 "class GSTLEARN_EXPORT Db: public AStringable, public ASerializable, public ICloneable"
    :return: Tuple (class name, list of inherited base classes)
    """
    # Modified regex to extract the class name and base classes, with or without braces on the same line
    pattern = r'class\s+\w+\s+(\w+)\s*[:\s]*([^{\n]+)?\s*(\{)?'

    match = re.match(pattern, line)

    if match:
        # Class name
        class_name = match.group(1)
        
        # Inherited base classes, separated by commas (if any)
        base_classes_raw = match.group(2)
        if base_classes_raw:
            base_classes_raw = base_classes_raw.split(',')
            # Remove the word 'public' from each base class
            base_classes = [base.replace('public', '').strip() for base in base_classes_raw]
        else:
            base_classes = []
        
        return class_name, base_classes
    else:
        return None, []

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
            if file.endswith('.hpp'):
                file_path = os.path.join(root, file)
                
                try:
                    with open(file_path, 'r', encoding='utf-8') as f:
                        for line in f:
                            # Extract class and base classes
                            class_name, base_classes = extract_class_and_bases(line)
                            
                            if class_name:
                                class_hierarchy[class_name] = base_classes
                                
                                # Check if the class inherits directly from AStringable
                                if 'AStringable' in base_classes:
                                    direct_inheritors.append(class_name)
                except (UnicodeDecodeError, FileNotFoundError):
                    continue

    # List of classes inheriting indirectly from AStringable
    all_inheritors = set(direct_inheritors)  # Start with direct inheritors

    # Traverse remaining classes to check for indirect inheritance
    newly_found = set(direct_inheritors)  # Classes already found

    # Iterate until no more new classes are found
    level = 2
    while newly_found:
        current_found = set()
        for class_name in class_hierarchy:
            if class_name not in all_inheritors:
                for base_class in class_hierarchy[class_name]:
                    if base_class in newly_found:
                        current_found.add(class_name)
                        break
        newly_found = current_found
        all_inheritors.update(newly_found)
        
        level += 1

    # Return the classes sorted alphabetically
    return all_inheritors

# Display all classes inheriting from AStringable (directly or indirectly) in alphabetical order

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


if __name__ == "__main__":
    print("--------------------------------------")
    output_txt_file = sys.argv[1]
    print(os.getcwd())
    print(output_txt_file)
    include_path = os.path.join("..", "..", "include")
    Astringable_classes = sorted(find_classes_inheriting_from_AStringable(include_path))
    with open(output_txt_file, "w", encoding='utf-8') as file:
        excluded = ["Node", "CovGradientFunctional", "GibbsUPropMono",
                    "Tapering", "GibbsUPropMultiMono", 
                    "CovGradientNumerical", "ElemNostat",
                    "GibbsMultiMono", "RuleShift", "GibbsUMultiMono",
                    "SpaceSN", "RuleShadow"]
        for class_name in Astringable_classes:
           if class_name not in excluded:
            file.write(generate_swig_extend_code(class_name))
