#!/usr/bin/env python3

import sys
import os
import utils
    
def replace_gstlrn_prefix(filename: str) -> None:
    with open(filename, 'r', encoding='utf-8') as f:
        content = f.read()

    # Ajouter la ligne de commentaire en tête
    content = "#File processed to clean namespace in default arguments\n" + content.replace('= gstlrn_', '= ')

    with open(filename, 'w', encoding='utf-8') as f:
        f.write(content)
#This script generates specific documentations when macros FORWARD_METHOD and FORWARD_METHOD_NON_CONST are used.


if __name__ == "__main__":
    folder = sys.argv[1]
    print(folder)
    replace_gstlrn_prefix(folder)
    
