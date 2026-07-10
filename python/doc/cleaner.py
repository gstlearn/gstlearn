import argparse
import os
import re
import xml.etree.ElementTree as ET


class DoxyGroup:
    def __init__(self, name, description, members):
        self.name = name
        self.description = description
        self.members = members


def clean_file(file_path):
    with open(file_path, "r", encoding="utf-8") as file:
        lines = file.readlines()

    text_to_insert = [
        "        if isinstance(self, AEnum):\n",
        '            strthis = self.getDescr() + " | " + str(self.getValue())\n',
        "            return strthis\n",
    ]
    lines[22:22] = text_to_insert

    cleaned_lines = []
    i = 0
    while i < len(lines):
        if "C++ includes:" in lines[i]:
            i += 1
            continue
        if i + 2 < len(lines) and 'r"""' in lines[i] and lines[i + 1].strip() == "":
            cleaned_lines.append(lines[i])
            cleaned_lines.append(lines[i + 1])
            i += 3
            continue
        cleaned_lines.append(lines[i])
        i += 1

    full_text = "".join(cleaned_lines)

    pattern = r"\$\s*([^$]*?)\s*\$"
    full_text = re.sub(pattern, r":math:`\1`", full_text)

    full_text = re.sub(r"\*\s*`(.*?)`", r"\1", full_text)

    with open(file_path, "w", encoding="utf-8") as file:
        file.write(full_text)


def parse_all_doxygen_groups(doxygen_xml_dir):
    all_groups = []
    index_path = os.path.join(doxygen_xml_dir, "index.xml")

    if not os.path.exists(index_path):
        raise FileNotFoundError()

    tree = ET.parse(index_path)
    root = tree.getroot()
    for compound in root.findall(".//compound[@kind='group']"):
        refid = compound.get("refid")
        group_xml_path = os.path.join(doxygen_xml_dir, f"{refid}.xml")

        if os.path.exists(group_xml_path):
            try:
                g_tree = ET.parse(group_xml_path)
                g_root = g_tree.getroot()
            except ET.ParseError:
                continue

            compounddef = g_root.find(".//compounddef[@kind='group']")
            if compounddef is None:
                continue

            title_node = compounddef.find("title")
            name = (
                title_node.text
                if title_node is not None
                else compounddef.find("compoundname").text
            )

            descriptions = []
            for para in compounddef.findall(".//briefdescription/para"):
                if para.text:
                    descriptions.append(para.text.strip())
            for para in compounddef.findall(".//detaileddescription/para"):
                if para.text:
                    descriptions.append(para.text.strip())

            description = "\n".join(descriptions)
            members = []
            for memberdef in compounddef.findall(".//memberdef"):
                name_node = memberdef.find("name")
                if name_node is not None and name_node.text:
                    members.append(name_node.text.strip())

            if description.strip():
                all_groups.append(DoxyGroup(name, description, members))

    return all_groups


def inject_grouped_docstrings(py_file_path, groups):
    member_to_group = {}
    group_members = {}
    group_definitions = {}

    for group in groups:
        group_definitions[group.name] = group.description.strip()
        group_members[group.name] = group.members
        for member in group.members:
            member_to_group[member] = group.name

    with open(py_file_path, "r", encoding="utf-8") as f:
        lines = f.read().splitlines(keepends=True)

    updated_lines = []
    processed_groups = set()
    skipped_functions = set()
    detected_classes = []

    i = 0
    while i < len(lines):
        line = lines[i]
        updated_lines.append(line)

        match = re.match(r"^(\s*)(def|class)\s+([a-zA-Z0-9_]+)\b", line)
        if match:
            type_ = match.group(2)
            indent = match.group(1)
            name = match.group(3)

            if type_ == "class" and not indent:
                detected_classes.append(name)

            if name in member_to_group:
                group_name = member_to_group[name]

                if group_name not in processed_groups:
                    description = group_definitions[group_name]
                    all_members = group_members[group_name]

                    doc_lines = []

                    for desc_line in description.split("\n"):
                        doc_lines.append(desc_line)

                    doc_lines.extend(
                        [
                            "",
                            ".. note:: Voir aussi les fonctions ci-dessous :",
                            "",
                            "   .. code-block:: none",
                            "",
                        ]
                    )

                    for m in all_members[1:]:
                        doc_lines.append(f"      {m}()")

                    doc_indent = indent + "    "
                    formatted_doc = (
                        f'{doc_indent}"""\n'
                        + "\n".join([f"{doc_indent}{l}".rstrip() for l in doc_lines])
                        + f'\n{doc_indent}"""\n'
                    )
                    updated_lines.append(formatted_doc)
                    processed_groups.add(group_name)
                else:
                    skipped_functions.add(name)
        i += 1

    with open(py_file_path, "w", encoding="utf-8") as f:
        f.writelines(updated_lines)

    return skipped_functions, detected_classes


def update_sphinx_conf(conf_py_path, skipped_functions):
    if not os.path.exists(conf_py_path):
        return

    with open(conf_py_path, "r", encoding="utf-8") as f:
        content = f.read()

    functions_list_str = str(list(skipped_functions))

    hook_code = f"""
def autodoc_skip_member_handler(app, what, name, obj, skip, options):
    short_name = name.split(".")[-1]

    if short_name == "thisown":
        return True

    if short_name.startswith("E_"):
        return True

    functions_to_skip = {functions_list_str}
    if name in functions_to_skip:
        return True
    return skip

def setup(app):
    app.connect('autodoc-skip-member', autodoc_skip_member_handler)
"""

    content = re.sub(
        r"def autodoc_skip_member_handler.*?app\.connect\('autodoc-skip-member', autodoc_skip_member_handler\)",
        "",
        content,
        flags=re.DOTALL,
    )
    content = content.strip() + "\n" + hook_code

    if "sphinx.ext.autosummary" not in content:
        content = content.replace(
            "extensions = [", "extensions = [\n    'sphinx.ext.autosummary',"
        )
    if "autosummary_generate" not in content:
        content += "\nautosummary_generate = True\n"

    with open(conf_py_path, "w", encoding="utf-8") as f:
        f.write(content)


def generate_autosummary_rst(output_rst_path, module_name, classes):
    rst_content = f"""Classes
=============

.. autosummary::
   :toctree: generated
   :template: class.rst

"""
    for c in sorted(classes):
        rst_content += f"   {module_name}.{c}\n"

    with open(output_rst_path, "w", encoding="utf-8") as f:
        f.write(rst_content)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("file_path")
    parser.add_argument("--xml-dir", required=True)
    parser.add_argument("--conf-py", required=True)
    parser.add_argument("--module-name", required=True)
    parser.add_argument("--output-rst", required=True)
    args = parser.parse_args()

    clean_file(args.file_path)
    groups = parse_all_doxygen_groups(args.xml_dir)
    skipped_funcs, detected_classes = inject_grouped_docstrings(args.file_path, groups)
    update_sphinx_conf(args.conf_py, skipped_funcs)
    generate_autosummary_rst(args.output_rst, args.module_name, detected_classes)
