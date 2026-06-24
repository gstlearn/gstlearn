#!/usr/bin/env python3
import re
from pathlib import Path


def md_to_rd(md_text):
    text = md_text.replace("\r\n", "\n")

    text = re.sub(
        r"(`{3,})[a-zA-Z0-9_-]*\s*\n(.*?)\n\1",
        r"\\preformatted{\n\2\n}",
        text,
        flags=re.S,
    )

    text = re.sub(
        r"('+)\s*\n(.*?)\n\1",
        r"\\preformatted{\n\2\n}",
        text,
        flags=re.S,
    )

    text = re.sub(r"[`']+\\preformatted\{", r"\\preformatted{", text)
    text = re.sub(r"\}[`']+", r"}", text)
    text = re.sub(r"^#+\s+(.*?)$", r"\n\\bold{\1}\n", text, flags=re.M)
    text = re.sub(r"\*\*(.*?)\*\*", r"\\bold{\1}", text)
    text = re.sub(r"\*(.*?)\*", r"\\emph{\1}", text)
    text = re.sub(r"`([^`\n]+)`", r"\\code{\1}", text)
    text = text.replace("\\rtools", "rtools")
    text = re.sub(r"\[(.*?)\]\((.*?)\)", r"\\href{\2}{\1}", text)
    text = re.sub(r"^[-\*]\s+(.*?)$", r"\\item \1", text, flags=re.M)
    text = re.sub(r"((?:\\item .*(?:\n|$))+)", r"\\itemize{\n\1}\n", text)
    text = re.sub(r"^\s*['`]+\s*$\n?", "", text, flags=re.M)

    return text


def generate_package_rd(readme_path, output_path):
    readme_content = Path(readme_path).read_text(encoding="utf-8")
    rd_details = md_to_rd(readme_content)

    rd_template = f"""\\name{{gstlearn}}
\\alias{{gstlearn}}
\\docType{{package}}
\\title{{gstlearn}}
\\details{{
{rd_details}
}}
"""

    output_file = Path(output_path)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    output_file.write_text(rd_template, encoding="utf-8")
    print(f"✅ Fichier .Rd nettoyé et valide généré dans : {output_path}")


if __name__ == "__main__":
    readme = "/home/mboisseau/a/gstlearn/r/README.md"
    output = "/home/mboisseau/a/gstlearn/build/r/Release/gstlearn/man/gstlearn.Rd"

    generate_package_rd(readme, output)
