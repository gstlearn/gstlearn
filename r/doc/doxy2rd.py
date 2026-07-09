#!/usr/bin/env python3
from enum import member
import os
from pathlib import Path
import argparse
import xml.etree.ElementTree as ET
import re
import csv


class Function:
    def __init__(
        self,
        name="",
        docType="",
        alias="",
        title="",
        description="",
        usage="",
        arguments=None,
        value="",
        format="",
        source="",
        references="",
        keyword="",
        note="",
        filename="",
    ):
        self.name = name
        self.filename = filename or name
        self.docType = docType
        self.alias = alias
        self.title = title
        self.description = description
        self.usage = usage
        self.arguments = arguments or []
        self.value = value
        self.format = format
        self.source = source
        self.references = references
        self.keyword = keyword
        self.note = note

    def to_rd(self):
        lines = []
        clean_name = (
            self.name.replace("gstlrn::", "")
            .replace("gstlrn_", "")
            .replace("[]", "_brackets")
        )
        if clean_name:
            lines.append(f"\\name{{{clean_name}}}")
            lines.append(f"\\alias{{{clean_name}}}")
        if self.docType:
            lines.append(f"\\docType{{{self.docType}}}")
        if self.title:
            lines.append(
                f"\\title{{{self.title.replace('gstlrn::', '').replace('gstlrn_', '').replace('::', '_')}}}"
            )
        if self.description:
            lines.append("\\description{\n" + self.description + "\n}")
        if self.usage:
            lines.append("\\usage{\n" + self.usage + "\n}")
        if self.arguments:
            lines.append("\\arguments{")
            for arg in self.arguments:
                lines.append(
                    f"  \\item{{{arg['name']}}}{{{arg.get('description', '')}}}"
                )
            lines.append("}")
        if self.value:
            lines.append("\\value{\n" + self.value + "\n}")
        if self.format:
            lines.append("\\format{\n" + self.format + "\n}")
        if self.note:
            lines.append("\\note{\n" + self.note + "\n}")
        if self.source:
            lines.append("\\source{\n" + self.source + "\n}")
        if self.references:
            lines.append("\\references{\n" + self.references + "\n}")
        if self.keyword:
            lines.append(f"\\keyword{{{self.keyword}}}")
        return "\n".join(lines)

    def write(self, output_dir):
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        safe_name = (
            self.filename.replace("gstlrn::", "")
            .replace("gstlrn_", "")
            .replace("/", "|")
            .replace("=", "_")
        )
        rd_file = output_dir / f"{safe_name}.Rd"
        with open(rd_file, "w", encoding="utf-8") as f:
            f.write(self.to_rd())
        return rd_file


class Class:
    def __init__(
        self,
        name="",
        title="",
        description="",
        methods=None,
        fields=None,
        groups=None,
        source="",
        references="",
        keyword="",
        note="",
    ):
        self.name = name
        self.title = title
        self.description = description
        self.methods = methods or []
        self.fields = fields or []
        self.groups = groups or []
        self.source = source
        self.references = references
        self.keyword = keyword
        self.note = note

    def to_rd(self):
        lines = []
        clean_name = (
            self.name.replace("gstlrn::", "")
            .replace("gstlrn_", "")
            .replace("[]", "_brackets")
        )
        if clean_name:
            lines.append(f"\\name{{{clean_name}}}")
            lines.append(f"\\alias{{{clean_name}}}")
        if self.title:
            lines.append(
                f"\\title{{{self.title.replace('gstlrn::', '').replace('gstlrn_', '')}}}"
            )
        if self.description:
            lines.append("\\description{\n" + self.description + "\n}")
        if self.fields:
            lines.append("\\format{")
            lines.append("  \\describe{")
            for f in self.fields:
                lines.append(f"    \\item{{{f['name']}}}{{{f.get('description', '')}}}")
            lines.append("  }")
            lines.append("}")

        if self.methods or self.groups:
            lines.append("\\details{")
            for g in self.groups:
                lines.append(f"\\strong{{{g['title']}}}")
                if g["description"]:
                    lines.append(f"\n{g['description']}\n")
                lines.append("  \\itemize{")

                for m_info in g["functions"]:
                    m_name = m_info["name"]
                    is_private = m_info["is_private"]
                    is_static = m_info["is_static"]
                    clean_m = (
                        m_name.replace("gstlrn::", "")
                        .replace("gstlrn_", "")
                        .replace("::", "_")
                    )
                    if is_private:
                        continue
                    if is_static:
                        lines.append(f"    \\item \\code{{{clean_name}_{clean_m}}}")
                    else:
                        lines.append(f"    \\item \\code{{{clean_m}}}")

                lines.append("  }\n")

            if self.methods:
                lines.append("Methods:")
                lines.append("  \\itemize{")

                class_prefix = f"{self.name.replace('gstlrn::', '').replace('gstlrn_', '').replace('::', '_')}_"

                for m in self.methods:
                    clean_m = (
                        m.replace("gstlrn::", "")
                        .replace("gstlrn_", "")
                        .replace("::", "_")
                    )

                    if clean_m.startswith(class_prefix):
                        display_name = clean_m[len(class_prefix) :]
                    else:
                        display_name = (
                            clean_m.split("_", 1)[-1] if "_" in clean_m else clean_m
                        )

                    display_name = re.sub(r"_\d+$", "", display_name)

                    lines.append(
                        f"    \\item \\link[={clean_m.replace('[]', '_brackets')}]{{{display_name}}}"
                    )
                lines.append("  }")
            lines.append("}")

        if self.note:
            lines.append("\\note{\n" + self.note + "\n}")
        if self.source:
            lines.append("\\source{" + self.source + "}")
        if self.references:
            lines.append("\\references{" + self.references + "}")
        if self.keyword:
            lines.append(f"\\keyword{{{self.keyword}}}")
        return "\n".join(lines)

    def write(self, output_dir):
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        safe_name = (
            self.name.replace("gstlrn::", "")
            .replace("gstlrn_", "")
            .replace("/", "|")
            .replace("=", "_")
        )
        rd_file = output_dir / f"{safe_name}.Rd"
        with open(rd_file, "w", encoding="utf-8") as f:
            f.write(self.to_rd())
        return rd_file


class Enumeration:
    def __init__(self, name, elements):
        self.name = name
        self.elements = elements

    def to_rd(self):
        text = f"\\name{{ {self.name} }}\n"
        text += f"\\alias{{ {self.name} }}\n"
        text += f"\\title{{ {self.name} }}\n"
        text += "\\details{{\n"
        text += "  \\tabular{lll}{\n"

        for idx, triplet in self.elements.items():
            text += f"    {triplet[0]} \\tab {triplet[1]} \\tab {triplet[2]} \\cr\n"

        text += "  }\n"
        text += "}}"
        return text

    def write(self, output_dir):
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        safe_name = (
            self.name.replace("gstlrn::", "")
            .replace("gstlrn_", "")
            .replace("/", "|")
            .replace("=", "_")
        )
        rd_file = output_dir / f"{safe_name}.Rd"
        with open(rd_file, "w", encoding="utf-8") as f:
            f.write(self.to_rd())
        return rd_file


class DoxygenParser:
    def __init__(self, xml_directory):
        self.xml_directory = Path(xml_directory)

    def _sanitize_rd_macros(self, text):
        text = re.sub(r"(?<!\\)%", r"\\%", text)
        text = text.replace("gstlrn::", "").replace("gstlrn_", "").replace("::", "_")

        valid_macros = {
            "name",
            "alias",
            "docType",
            "title",
            "description",
            "usage",
            "arguments",
            "item",
            "value",
            "format",
            "note",
            "source",
            "references",
            "keyword",
            "details",
            "itemize",
            "describe",
            "code",
            "link",
            "strong",
            "eqn",
            "deqn",
            "enumerate",
        }

        def replacer(match):
            macro_name = match.group(1)
            if macro_name not in valid_macros:
                return "|" + macro_name
            return match.group(0)

        return re.sub(r"\\([a-zA-Z]+)", replacer, text)

    def _text(self, node, skip_params=False, skip_returns=False, skip_remarks=False):
        if node is None:
            return ""
        parts = []
        state = {"in_listitem": False}

        def get_inner_text(elem):
            res = [elem.text or ""]
            for c in elem:
                res.append(get_inner_text(c))
                res.append(c.tail or "")
            return "".join(res)

        def walk(elem):
            if (
                skip_remarks
                and elem.tag == "simplesect"
                and elem.attrib.get("kind") in ("remark", "note")
            ):
                return
            if (
                skip_returns
                and elem.tag == "simplesect"
                and elem.attrib.get("kind") == "return"
            ):
                return
            if (
                skip_params
                and elem.tag == "parameterlist"
                and elem.attrib.get("kind") == "param"
            ):
                return

            if elem.tag == "itemizedlist":
                parts.append("\n\\itemize{\n")

            if elem.tag == "orderedlist":
                parts.append("\n\\enumerate{\n")

            if elem.tag == "listitem":
                parts.append("  \\item ")
                state["in_listitem"] = True

            if elem.tag == "ref":
                text_content = get_inner_text(elem).strip()
                clean_name = (
                    text_content.replace("gstlrn::", "")
                    .replace("gstlrn_", "")
                    .replace("()", "")
                    .replace("::", "_")
                    .strip()
                )
                clean_text = text_content.replace("gstlrn::", "").replace("gstlrn_", "")
                parts.append(f"\\link[={clean_name}]{{{clean_text}}}")
                return

            if elem.tag == "formula":
                form_text = (elem.text or "").strip()
                if form_text.startswith("\\[") and form_text.endswith("\\]"):
                    eq = form_text[2:-2].strip()
                    parts.append(f"\n\\deqn{{{eq}}}{{{eq}}}\n")
                elif form_text.startswith("$") and form_text.endswith("$"):
                    eq = form_text[1:-1].strip()
                    parts.append(f"\\eqn{{{eq}}}{{{eq}}}")
                else:
                    parts.append(f"\\eqn{{{form_text}}}{{{form_text}}}")
                return

            if elem.tag == "computeroutput" and elem.text and "<br/>" in elem.text:
                parts.append(elem.text.replace("<br/>", "\n\n"))
                return

            if elem.tag in ("para", "simplesect"):
                if not state["in_listitem"] and parts:
                    parts.append("\n\n")

            if elem.text:
                clean_text = re.sub(r"[\r\n\t]+", " ", elem.text)
                parts.append(clean_text)

            for child in elem:
                if child.tag == "linebreak" or child.tag == "br":
                    parts.append("\n\n")

                walk(child)

                if child.tail:
                    clean_tail = re.sub(r"[\r\n\t]+", " ", child.tail)
                    parts.append(clean_tail)

            if elem.tag == "listitem":
                parts.append("\n")
                state["in_listitem"] = False

            if elem.tag == "orderedlist":
                parts.append("}\n")

            if elem.tag == "itemizedlist":
                parts.append("}\n")

        walk(node)

        full_text = "".join(parts)
        full_text = (
            full_text.replace("<br/>", "\n\n")
            .replace("<br >", "\n\n")
            .replace("<br>", "\n\n")
        )
        full_text = re.sub(r"[ \t]+", " ", full_text)
        full_text = re.sub(r"\n {1,}", "\n", full_text)
        full_text = re.sub(r"\n{3,}", "\n\n", full_text)

        return self._sanitize_rd_macros(full_text.strip())

    def _extract_returns(self, node):
        if node is None:
            return ""
        returns = []
        for simplesect in node.findall(".//simplesect"):
            if simplesect.attrib.get("kind") == "return":
                returns.append(self._text(simplesect))
        return self._format_links("\n\n".join(returns).strip())

    def _extract_remarks(self, node):
        if node is None:
            return ""
        remarks = []
        for simplesect in node.findall(".//simplesect"):
            if simplesect.attrib.get("kind") in ("remark", "note"):
                remarks.append(self._text(simplesect))
        return self._format_links("\n\n".join(remarks).strip())

    def _format_links(self, text):
        pattern = r"(?:\\|@)link\s+(\S+)\s+(.*?)\s*(?:\\|@)endlink"

        def link_replacer(match):
            target = (
                match.group(1)
                .replace("gstlrn::", "")
                .replace("gstlrn_", "")
                .replace("::", "_")
            )
            display = match.group(2).replace("gstlrn::", "").replace("gstlrn_", "")
            return f"\\link[={target}]{{{display}}}"

        text = re.sub(pattern, link_replacer, text)
        return self._sanitize_rd_macros(text)

    def _get_arguments(self, memberdef):
        arguments = []
        for param in memberdef.findall("param"):
            pname = param.findtext("declname", "")
            arguments.append({"name": pname, "description": ""})

        detailed = memberdef.find("detaileddescription")
        if detailed is not None:
            for plist in detailed.findall(".//parameterlist"):
                if plist.attrib.get("kind") != "param":
                    continue
                for item in plist.findall("parameteritem"):
                    names = [self._text(x) for x in item.findall(".//parametername")]
                    desc = self._text(item.find("parameterdescription"))
                    desc = self._format_links(desc)
                    for arg in arguments:
                        if arg["name"] in names:
                            arg["description"] = desc
        return arguments

    def parse(self):
        functions = []
        classes = []

        doxygen_groups = {}
        for group_file in self.xml_directory.glob("group__*.xml"):
            try:
                g_tree = ET.parse(group_file)
                g_root = g_tree.getroot()
                for compdef in g_root.findall(".//compounddef"):
                    g_id = compdef.attrib.get("id", "")
                    g_title = compdef.findtext("title", g_id)

                    g_brief = self._text(
                        compdef.find("briefdescription"),
                        skip_params=True,
                        skip_returns=True,
                        skip_remarks=True,
                    )
                    g_detailed = self._text(
                        compdef.find("detaileddescription"),
                        skip_params=True,
                        skip_returns=True,
                        skip_remarks=True,
                    )
                    g_desc = self._format_links(
                        "\n\n".join(x for x in [g_brief, g_detailed] if x)
                    )

                    g_funcs = []
                    for m_node in compdef.findall(".//memberdef[@kind='function']"):
                        g_m_name = m_node.findtext("name", "")
                        if any(
                            word in g_m_name
                            for word in ["operator", "~", "final_scalar_type"]
                        ):
                            continue

                        g_funcs.append(
                            {
                                "name": m_node.findtext("name", ""),
                                "is_private": m_node.attrib.get("prot")
                                in ("private", "protected"),
                                "is_static": m_node.attrib.get("static") == "yes",
                            }
                        )

                    doxygen_groups[g_id] = {
                        "title": g_title,
                        "description": g_desc,
                        "functions": g_funcs,
                    }
            except Exception:
                continue

        for xml_file in self.xml_directory.glob("*.xml"):
            if xml_file.name.startswith("group__"):
                continue
            try:
                tree = ET.parse(xml_file)
                root = tree.getroot()
            except Exception:
                continue

            for compounddef in root.findall(".//compounddef"):
                kind = compounddef.attrib.get("kind", "")

                if kind not in ("class", "struct", "file", "namespace"):
                    continue

                c_name = compounddef.findtext("compoundname", "")
                brief_node = compounddef.find("briefdescription")
                detailed_node = compounddef.find("detaileddescription")

                brief = self._text(
                    brief_node, skip_params=True, skip_returns=True, skip_remarks=True
                )
                detailed = self._text(
                    detailed_node,
                    skip_params=True,
                    skip_returns=True,
                    skip_remarks=True,
                )
                desc = self._format_links(
                    "\n\n".join(x for x in [brief, detailed] if x)
                )

                c_remark = "\n\n".join(
                    x
                    for x in [
                        self._extract_remarks(brief_node),
                        self._extract_remarks(detailed_node),
                    ]
                    if x
                )

                methods = []
                fields = []
                class_bound_groups = {}
                method_counts = {}

                for member_node in compounddef.findall(".//member"):
                    refid = member_node.attrib.get("refid", "")
                    if (
                        refid.startswith("group__")
                        and member_node.attrib.get("kind") == "function"
                    ):
                        group_key = refid.split("_1g")[0] if "_1g" in refid else refid
                        if group_key not in doxygen_groups:
                            possible_keys = [
                                k for k in doxygen_groups if refid.startswith(k)
                            ]
                            group_key = possible_keys[0] if possible_keys else None

                        if group_key and group_key in doxygen_groups:
                            if group_key not in class_bound_groups:
                                class_bound_groups[group_key] = doxygen_groups[
                                    group_key
                                ]

                for section in compounddef.findall("sectiondef"):
                    for member in section.findall("memberdef"):
                        if member.attrib.get("prot") in ("private", "protected"):
                            continue
                        m_name = member.findtext("name", "")
                        m_kind = member.attrib.get("kind", "")

                        if m_kind == "function":
                            skip_words = ["operator", "~", "final_scalar_type"]
                            is_in_any_group = any(
                                any(f["name"] == m_name for f in g["functions"])
                                for g in class_bound_groups.values()
                            )
                            if is_in_any_group:
                                continue

                            if kind in ("file", "namespace"):
                                f_name = m_name
                            else:
                                f_name = f"{c_name}_{m_name}"

                            if any(word in m_name for word in skip_words):
                                continue

                            if m_name not in method_counts:
                                method_counts[m_name] = 1
                                f_filename = f_name
                            else:
                                method_counts[m_name] += 1
                                f_filename = f"{f_name}_{method_counts[m_name]}"

                            if kind in ("class", "struct"):
                                methods.append(f_filename)

                            f_detailed_node = member.find("detaileddescription")
                            f_desc = self._format_links(
                                self._text(
                                    f_detailed_node,
                                    skip_params=True,
                                    skip_returns=True,
                                    skip_remarks=True,
                                )
                            )

                            title = m_name
                            usage_f_name = m_name
                            if member.attrib.get("static") == "yes" and kind in (
                                "class",
                                "struct",
                            ):
                                title = f_name
                                usage_f_name = f_name

                            functions.append(
                                Function(
                                    name=f_name,
                                    title=title,
                                    description=f_desc,
                                    usage=f"{usage_f_name}{member.findtext('argsstring', '')}",
                                    arguments=self._get_arguments(member),
                                    value=self._extract_returns(f_detailed_node),
                                    note=self._extract_remarks(f_detailed_node),
                                    filename=f_filename,
                                )
                            )
                        elif m_kind in ("variable", "property") and kind in (
                            "class",
                            "struct",
                        ):
                            f_field_desc = self._format_links(
                                self._text(member.find("briefdescription"))
                            )
                            fields.append({"name": m_name, "description": f_field_desc})

                if kind in ("class", "struct"):
                    if any(word in c_name for word in skip_words):
                        continue

                    classes.append(
                        Class(
                            name=c_name,
                            title=c_name,
                            description=desc,
                            methods=methods,
                            fields=fields,
                            groups=list(class_bound_groups.values()),
                            note=c_remark,
                        )
                    )

        return functions, classes

    def enum_parser(self, xml_file):
        tree = ET.parse(xml_file)
        root = tree.getroot()

        initializer_node = root.find(".//memberdef[@kind='define']/initializer")

        if initializer_node is not None:
            raw_text = "".join(initializer_node.itertext())
            clean_text = raw_text.replace("\\\n", "").replace("\\", "").strip()

            reader = csv.reader([clean_text], skipinitialspace=True)
            elements = next(reader)

            enum_name = elements[0]

            enum_dict = {}
            counter = 1

            for i in range(2, len(elements), 3):
                if i + 2 < len(elements):
                    key = elements[i]
                    value = elements[i + 1]
                    desc = elements[i + 2].strip()

                    enum_dict[counter] = [key, value, desc]
                    counter += 1

            return Enumeration(enum_name, enum_dict)
        else:
            return None


if __name__ == "__main__":
    parser_args = argparse.ArgumentParser()
    parser_args.add_argument("INPUT_DIR")
    parser_args.add_argument("OUTPUT_DIR")
    args = parser_args.parse_args()

    print(f"Parsing Doxygen XML files from: {args.INPUT_DIR}")
    parser = DoxygenParser(args.INPUT_DIR)
    funcs, classes = parser.parse()

    print(f"Writing output files to: {args.OUTPUT_DIR}")
    for f in funcs:
        f.write(args.OUTPUT_DIR)
    for c in classes:
        c.write(args.OUTPUT_DIR)

    print(f"Parsing Enumeration XML files from: {args.INPUT_DIR}")
    for file in os.listdir(args.INPUT_DIR):
        if file.endswith(".xml") and file.startswith("E"):
            xml_file_path = os.path.join(args.INPUT_DIR, file)
            try:
                enum_obj = parser.enum_parser(xml_file_path)
                if enum_obj is not None:
                    enum_obj.write(args.OUTPUT_DIR)
            except Exception as e:
                print(f"Error parsing {file}: {e}")
                continue
