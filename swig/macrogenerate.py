import re


def extract_macro_calls(header_file = "../include/Model/ModelGeneric.hpp", macro_name = "FORWARD_METHOD"):
    """
    Extrait les appels à une macro spécifique dans un fichier d'en-tête.
    """
    macro_calls = []
    macro_pattern =  rf"^(?!#define\s+{macro_name})\s*{macro_name}\(([^,]+),\s*([^)]*)\)"
    with open(header_file, 'r') as file:
        for line in file:
            match = re.search(macro_pattern, line)
            if match:
                print(line)
                body, func_name = match.groups()
                macro_calls.append((func_name, body.strip()))
    return macro_calls

def generate_python_code(macro_calls, output_file, classname):
    """
    Génère du code Python pour les macros extraites.
    """
    
    with open(output_file, 'w') as file:
        file.write('%pythoncode %{\n')
        file.write(f"import gstlearn as gl\n")
        for func_name, body in macro_calls:
            # Conversion simple : remplacer le corps C++ par une chaîne Python équivalente
            python_body = body.replace('std::cout', 'print').replace(';', '')
            file.write("\n")
            file.write(f"def {func_name}(self,*args,**kwargs):\n")
            file.write(f"    return self.{python_body}().{func_name}(*args,**kwargs)\n")
            file.write("\n")
            file.write(f"setattr(gl.{classname},'{func_name}',{func_name})\n")
        file.write('%}\n')


def build_macro_python(macro_name, filename,output_file = "../../python/generated_python.i"):
  header_file = "../../include/"+ filename
  macro_calls = extract_macro_calls(header_file, macro_name)
  generate_python_code(macro_calls, output_file, filename.split("/")[0])
  print(f"Code Python généré dans {output_file}")

def build_macro_r(macro_name, filename,output_file = "../../python/generated_python.i"):
  header_file = "../../include/"+ filename
  macro_calls = extract_macro_calls(header_file, macro_name)
  generate_r_code(macro_calls, output_file, filename.split("/")[0])
  print(f"Code R généré dans {output_file}")

header_file = "Model/ModelGeneric.hpp"  # Fichier d'en-tête à analyser
macro_name = "FORWARD_METHOD"  # Nom de la macro à analyser
build_macro_python(macro_name,header_file)




