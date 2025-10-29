"""
create Latex PDF from jupytext in notebooks/
"""

import shutil
import subprocess
import sys
from pathlib import Path

import yaml


def run(cmd, cwd=None):
    """Run a shell command with error checking."""
    print(f"\n▶ Running: {cmd}")
    result = subprocess.run(cmd, shell=True, cwd=cwd, check=False)
    if result.returncode != 0:
        sys.exit(f"❌ Command failed: {cmd}")


def convert_jupytext(arg):
    """Convert jupytext to md file for using it in jupyter-book"""
    # Source and target directories
    py_file = str(arg) + ".py"
    output_dir = Path("book")

    # Make sure output directory exists
    output_dir.mkdir(parents=True, exist_ok=True)

    base = arg  # filename without extension
    output_file = output_dir / f"{base}.md"
    # Run the jupytext command
    subprocess.run(
        [
            "uv",
            "run",
            "jupytext",
            "--to",
            "myst",
            "notebooks/" + py_file,
            "--output",
            str(output_file),
        ],
        check=True,
    )

    print(f"Converted {py_file} -> {output_file}")


def update_toc(toc_path, arg):
    """Update the file entry in _toc.yml with the given argument."""
    toc_path = Path(toc_path)
    toc = yaml.safe_load(toc_path.read_text(encoding="utf-8"))

    # Traverse sections/chapters and replace the file name
    def replace_file(node):
        if isinstance(node, list):
            for item in node:
                replace_file(item)
        elif isinstance(node, dict):
            if "file" in node:
                node["file"] = arg
            if "sections" in node:
                replace_file(node["sections"])
            if "chapters" in node:
                replace_file(node["chapters"])

    replace_file(toc.get("sections") or toc.get("chapters") or [])
    toc_path.write_text(yaml.dump(toc, sort_keys=False), encoding="utf-8")
    print(f"✅ Updated TOC file to use: {arg}")


def update_config(config_path, title):
    """Update the 'title' field in _config.yml."""
    config_path = Path(config_path)
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))

    # Update title (and ensure key exists)
    config["title"] = title

    config_path.write_text(yaml.dump(config, sort_keys=False), encoding="utf-8")
    print(f"✅ Updated book title in config to: {title}")


def build_pdf(arg, title):
    """Build PDF files"""
    project_root = Path(__file__).parent.parent
    toc_path = project_root / "book" / "_toc.yml"
    config_path = project_root / "book" / "_config.yml"
    latex_dir = project_root / "book" / "_build" / "latex"
    output_pdf = project_root / "book" / "deliverables" / f"{arg}.pdf"

    # Step 1: update md documents
    convert_jupytext(arg)

    # Step 2: overwrite toc
    update_toc(toc_path, arg)

    # Step 3: update title in config
    update_config(config_path, title)

    # Step 4: build latex docs
    run("uv run jupyter-book build book/ --builder latex", cwd=project_root)

    # Step 5: build latex pdf
    tex_file = latex_dir / "projectnamenotset.tex"
    if not tex_file.exists():
        sys.exit(f"❌ TeX file not found: {tex_file}")

    with open(tex_file) as fh:
        tex_raw = fh.read()

    tex_raw = tex_raw.replace(
        r"\end{document}",
        "\\printbibliography[heading=bibintoc,title={Whole bibliography}]\n\\end{document}",
    )
    with open(tex_file, "w") as fh:
        fh.write(tex_raw)

    run(f"xelatex {tex_file.name}", cwd=latex_dir)
    shutil.copy(project_root / "book" / "references.bib", latex_dir / "references.bib")
    run(f"biber {tex_file.name.replace('.tex', '.bcf')}", cwd=latex_dir)
    # second pass for updating cross-references
    run(f"xelatex {tex_file.name}", cwd=latex_dir)
    run(f"xelatex {tex_file.name}", cwd=latex_dir)

    # Step 6: copy result
    built_pdf = latex_dir / "projectnamenotset.pdf"
    if not built_pdf.exists():
        sys.exit("❌ PDF was not generated.")
    output_pdf.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy(built_pdf, output_pdf)
    print(f"✅ Saved final PDF → {output_pdf}")


if __name__ == "__main__":
    if len(sys.argv) != 3:  # noqa: PLR2004
        sys.exit("Usage: python create_latex_pdf.py <arg> <title>")
    arg = sys.argv[1]
    title = sys.argv[2]
    build_pdf(arg, title)
