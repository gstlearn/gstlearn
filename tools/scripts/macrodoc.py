#!/usr/bin/env python3

import sys
import os
import utils
    
    
#This script generates specific documentations when macros FORWARD_METHOD and FORWARD_METHOD_NON_CONST are used.

if __name__ == "__main__":
   folder  = sys.argv[1]
   for macro_name in ["FORWARD_METHOD","FORWARD_METHOD_NON_CONST"]:
    tree = os.walk(os.path.join(os.pardir, "include"))
    for dir, root, files in tree:
        path = os.path.split(dir)[1]
        if path in ["include", "Core"]:
            continue
        for file in files:
            if file.split(".")[-1] != "hpp":
                continue
            header_file = os.path.join(os.pardir,"include",path, file)  # Utilisation de os.path.join
            
            founds = utils.extract_macro_calls(header_file, macro_name)
            if len(founds) > 0 :
                output_file = os.path.join(folder, file)
                with open(output_file, "w") as ofile:
                    ofile.write("class " + file.split(".")[0] + "{\n")
                    ofile.write("public:\n")
                    for found in founds:
                        method_get = found[2]
                        classname = utils.find_function_return_type_in_file(method_get, header_file)
                        groupfolder = utils.find_include_folder_in_file(classname, header_file)
                        method_name = found[1]
                        header_path = os.path.join(os.pardir,"include",groupfolder,classname + ".hpp")
                        signature = utils.find_method_signature_in_header(method_name, header_path)

                        
                        
                        if signature is None:
                            signature = utils.find_method_in_class_or_base(classname, method_name, os.path.join(os.pardir,"include"))
                            if signature is None:
                                print(f"⚠️ Méthode non trouvée ({file}): {method_name} dans {header_path}")
                                continue
                         
           
                        ofile.write("/**\n")
                        ofile.write("* @brief Call the method \\ref " + classname + "::" + method_name +" of the object \\ref " + classname + ".\n")
                        ofile.write("*/\n")
                        ofile.write(signature +"\n")
                        
                        
                        
                    ofile.write("};")
       
