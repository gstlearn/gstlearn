#!/usr/bin/env python3

import sys
import os
import utils
    
    
#This script generates specific documentations when macros FORWARD_METHOD is used.

if __name__ == "__main__":
   folder  = sys.argv[1]
   for macro_name in ["FORWARD_METHOD","FORWARD_METHOD_NON_CONST"]:
    tree = os.walk(os.path.join(os.pardir, os.pardir, "include"))
    for dir, root, files in tree:
        path = os.path.split(dir)[1]
        if path in ["include", "Core"]:
            continue
        for file in files:
            if file.split(".")[-1] != "hpp":
                continue
            header_file = os.path.join(path, file)  # Utilisation de os.path.join
            print(extrac_macro_calls(header_file, macro_name))
   print("coucou from " + folder)
