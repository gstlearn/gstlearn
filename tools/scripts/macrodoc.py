#!/usr/bin/env python3

import sys
import os
import utils
    
    
#This script generates specific documentations when macros FORWARD_METHOD and FORWARD_METHOD_NON_CONST are used.


if __name__ == "__main__":
    folder = sys.argv[1]
    signatures_by_file = {}

    for macro_name in ["FORWARD_METHOD", "FORWARD_METHOD_NON_CONST"]:
        for dir, _, files in os.walk(os.path.join(os.pardir, "include")):
            path = os.path.split(dir)[1]
            if path in ["include", "Core"]:
                continue

            for file in files:
                if not file.endswith(".hpp"):
                    continue

                header_file = os.path.join(os.pardir, "include", path, file)
                founds = utils.extract_macro_calls(header_file, macro_name)
                if not founds:
                    continue

                class_name = file.split(".")[0]
                output_file = os.path.join(folder, file)

                if file not in signatures_by_file:
                    signatures_by_file[file] = {
                        "class_name": class_name,
                        "signatures": []
                    }

                for found in founds:
                    method_name = found[1]
                    method_get = found[2]

                    classname = utils.find_function_return_type_in_file(method_get, header_file)
                    method_file = utils.find_header_for_class(classname, os.path.join(os.pardir, "include"))
                    if not method_file:
                        continue

                    signature = utils.find_method_signature_in_header(method_name, method_file)
                    if signature is None:
                        signature = utils.find_method_in_class_or_bases(classname, method_name, os.path.join(os.pardir, "include"))
                        if signature is None:
                            print(f"⚠️ Méthode non trouvée ({file}): {method_name}")
                            continue

                    doc = f"/**\n * @brief Call the method \\ref {classname}::{method_name} of the object \\ref {classname}.\n */"
                    signatures_by_file[file]["signatures"].append(doc + "\n" + signature)

    # Génération finale
    for file, content in signatures_by_file.items():
        output_path = os.path.join(folder, file)
        with open(output_path, "w") as f:
            f.write(f"class {content['class_name']} {{\npublic:\n")
            for sig in content["signatures"]:
                f.write(sig + "\n\n")
            f.write("};\n")
