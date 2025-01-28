#!/usr/bin/env python3
import re
import sys
import os

def extract_macro_calls(header_file , macro_name):
    classname = header_file.split("/")[3].split(".")[0]
    macro_calls = []
    macro_pattern = macro_pattern = rf"^(?!#define\s+{macro_name})\s*{macro_name}\((.*?)\)"
    with open(header_file, 'r') as file:
        for line in file:
            match = re.findall(macro_pattern, line, re.MULTILINE)
            if match:
                contents = match[0].split(",")
                if len(contents) == 2:
                    body, func_name = contents
                    arg = ""
                else :
                    if len(contents) == 3:
                        body, func_name,arg = contents
                    else:
                        print("error.")  
                macro_calls+= [(classname, func_name.strip(), body.strip(),arg)]
    return macro_calls
def generate_python_code(macro_calls, output_file,first = True):
    """
    Génère du code Python pour les macros extraites.
    """
    if len(macro_calls) == 0:
        return first
    if first:
        opt = "w"
    else:
        opt = "a"
    with open(output_file, opt) as file:
        file.write('%pythoncode %{\n')
        file.write(f"import gstlearn as gl\n")
        for classname, func_name, body, arg in macro_calls:
            # Conversion simple : remplacer le corps C++ par une chaîne Python équivalente
            python_body = body.replace('std::cout', 'print').replace(';', '')
            file.write("\n")
            file.write(f"def {func_name}(self,*args,**kwargs):\n")
            file.write(f"    return self.{python_body}().{func_name}(*args,**kwargs)\n")
            file.write("\n")
            file.write(f"setattr(gl.{classname},'{func_name}',{func_name})\n")
        file.write('%}\n')
    return False

def generate_r_code(macro_calls, output_file,first = True):
    
    if len(macro_calls) == 0:
        return first
    if first:
        opt = "w"
    else:
        opt = "a"
    with open(output_file, opt) as file:
        file.write("%insert(s)%{\n")

        for classname, func_name, body, arg in macro_calls:
            python_body = body.replace('std::cout', 'print').replace(';', '')
            file.write(f"f = function(self,...)\n")
            file.write(f"{{\n")
            file.write(f"   return({classname}_{python_body}(self)${func_name}(...))\n")
            file.write(f"}}\n")
            file.write(f"\n")
            file.write(f"assign('{classname}_evalCovMatSym', f , envir = asNamespace('gstlearn'))\n")
        
        file.write("%}\n")
    return False
     
def build_macro_python(macro_name, filename,output_file, first = True):
  header_file = "../../include/"+ filename
  output_file = "../../python/" + output_file
  macro_calls = extract_macro_calls(header_file, macro_name)
  first = generate_python_code(macro_calls, output_file, first)
  print(f"Python additional wrappers generated in {output_file}")
  return first

def build_macro_r(macro_name, filename,output_file,first = True):
  header_file = "../../include/"+ filename
  output_file = "../../r/" + output_file
  macro_calls = extract_macro_calls(header_file, macro_name)
  first = generate_r_code(macro_calls, output_file, first)
  print(f"R additional wrappers generated in {output_file}")
  return first

macro_name = "FORWARD_METHOD"  # Nom de la macro à analyser
header_file = "Model/ModelGeneric.hpp"  # Fichier d'en-tête à analyser
output_file = "generated.i"  # Fichier de sortie
includedir = "../../include/"

tree = os.walk(os.path.join(os.pardir,os.pardir,"include"))
first = True
for dir,root,files in tree:
    path = os.path.split(dir)[1]
    if path in ["include","Core"]:
        continue
    for file in files:
        if file.split(".")[-1] != "hpp":
            continue
        header_file = os.path.join(path,file)
        if sys.argv[1] == "python":
            first = build_macro_python(macro_name,header_file,"generated_python.i",first)
        if sys.argv[1] == "r":
            first = build_macro_r(macro_name,header_file,"generated_r.i",first)



