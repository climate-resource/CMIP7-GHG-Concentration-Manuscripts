#!/usr/bin/env python3
"""Exercise a real exported manuscript's assisted return path on disposable files."""

import argparse
import copy
import importlib.util
import json
import re
import shutil
import subprocess
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path
from types import SimpleNamespace
from zipfile import ZipFile

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    "review", ROOT / "scripts/manuscript-word-review.py"
)
review = importlib.util.module_from_spec(spec)
spec.loader.exec_module(review)
W, NS = review.W, review.NS


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manuscript", choices=("historical", "scenario"), required=True
    )
    parser.add_argument(
        "--compile",
        action="store_true",
        help="Also recompile the edited disposable LaTeX and run its value checks",
    )
    args = parser.parse_args()
    name = args.manuscript + "-ghg-forcing-for-cmip7"
    manifests = sorted((ROOT / "review-baselines" / name).glob("*/manifest.json"))
    if not manifests:
        raise SystemExit("Export this manuscript first.")
    manifest_path = manifests[-1]
    manifest = json.loads(manifest_path.read_text())
    intro_path = f"manuscripts/{name}/introduction.tex"
    intro = manifest["snapshots"][intro_path]["text"]
    with tempfile.TemporaryDirectory(
        prefix="mm-review-test-", dir=ROOT / "build"
    ) as temporary:
        tmp = Path(temporary)
        for path, data in manifest["snapshots"].items():
            target = tmp / path
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_text(data["text"])
        baseline_dir = tmp / "baseline"
        baseline_dir.mkdir()
        shutil.copy2(manifest_path, baseline_dir / "manifest.json")
        shutil.copy2(
            manifest_path.parent / "baseline.docx", baseline_dir / "baseline.docx"
        )
        with ZipFile(baseline_dir / "baseline.docx") as archive:
            entries = {n: archive.read(n) for n in archive.namelist()}
        document = ET.fromstring(entries["word/document.xml"])
        paragraphs = document.findall(".//w:body//w:p", NS)
        paragraph_info = {
            p["bookmark"]: p for p in manifest["paragraphs"] if p["bookmark"]
        }
        candidates = []
        for p in paragraphs:
            names = [b.get(f"{{{W}}}name") for b in p.findall(".//w:bookmarkStart", NS)]
            hint = next((paragraph_info[n] for n in names if n in paragraph_info), {})
            if hint.get("source") == intro_path:
                candidates.append(p)
        tracked = next(
            p
            for p in candidates
            if any("seventh" in (t.text or "") for t in p.findall(".//w:t", NS))
        )
        original_math = [ET.tostring(m) for m in document.findall(".//m:oMath", NS)]
        text_node = next(
            t for t in tracked.findall(".//w:t", NS) if "seventh" in (t.text or "")
        )
        original_run = next(
            r for r in tracked.findall("w:r", NS) if text_node in list(r)
        )
        position = list(tracked).index(original_run)
        deletion = ET.Element(
            f"{{{W}}}del", {f"{{{W}}}id": "70001", f"{{{W}}}author": "MM Test"}
        )
        deleted_run = copy.deepcopy(original_run)
        for t in deleted_run.findall("w:t", NS):
            t.tag = f"{{{W}}}delText"
        deletion.append(deleted_run)
        insertion = ET.Element(
            f"{{{W}}}ins", {f"{{{W}}}id": "70002", f"{{{W}}}author": "MM Test"}
        )
        inserted_run = copy.deepcopy(original_run)
        inserted_run.find("w:t", NS).text = text_node.text.replace("seventh", "7th", 1)
        insertion.append(inserted_run)
        tracked.remove(original_run)
        tracked.insert(position, deletion)
        tracked.insert(position + 1, insertion)
        ordinary = None
        for p in candidates:
            if p is tracked:
                continue
            for t in p.findall(".//w:t", NS):
                phrase = t.text or ""
                if (
                    re.fullmatch(r"[A-Za-z][A-Za-z ,'-]{14,}", phrase)
                    and intro.count(phrase) == 1
                ):
                    ordinary = p
                    old_phrase = phrase
                    t.text = phrase + " review verification"
                    break
            if ordinary is not None:
                break
        assert ordinary is not None, "No safe ordinary-edit test anchor found"
        sentence = "This paragraph verifies the Word review transfer workflow."
        new_paragraph = ET.Element(f"{{{W}}}p")
        ET.SubElement(
            ET.SubElement(new_paragraph, f"{{{W}}}r"), f"{{{W}}}t"
        ).text = sentence
        body = document.find("w:body", NS)
        body.insert(list(body).index(ordinary) + 1, new_paragraph)
        start = ET.Element(f"{{{W}}}commentRangeStart", {f"{{{W}}}id": "70003"})
        end = ET.Element(f"{{{W}}}commentRangeEnd", {f"{{{W}}}id": "70003"})
        tracked.insert(0, start)
        tracked.append(end)
        reference_run = ET.SubElement(tracked, f"{{{W}}}r")
        ET.SubElement(
            reference_run, f"{{{W}}}commentReference", {f"{{{W}}}id": "70003"}
        )
        comment_root = ET.Element(f"{{{W}}}comments")
        comment = ET.SubElement(
            comment_root,
            f"{{{W}}}comment",
            {f"{{{W}}}id": "70003", f"{{{W}}}author": "MM Test"},
        )
        ET.SubElement(
            ET.SubElement(ET.SubElement(comment, f"{{{W}}}p"), f"{{{W}}}r"), f"{{{W}}}t"
        ).text = "Please review this citation context."
        entries["word/comments.xml"] = ET.tostring(
            comment_root, encoding="utf-8", xml_declaration=True
        )
        relns = "http://schemas.openxmlformats.org/package/2006/relationships"
        rels = ET.fromstring(entries["word/_rels/document.xml.rels"])
        ET.SubElement(
            rels,
            f"{{{relns}}}Relationship",
            {
                "Id": "rIdMMComment",
                "Type": "http://schemas.openxmlformats.org/officeDocument/2006/relationships/comments",
                "Target": "comments.xml",
            },
        )
        entries["word/_rels/document.xml.rels"] = ET.tostring(rels)
        ct = ET.fromstring(entries["[Content_Types].xml"])
        ET.SubElement(
            ct,
            "{http://schemas.openxmlformats.org/package/2006/content-types}Override",
            {
                "PartName": "/word/comments.xml",
                "ContentType": "application/vnd.openxmlformats-officedocument."
                "wordprocessingml.comments+xml",
            },
        )
        entries["[Content_Types].xml"] = ET.tostring(ct)
        entries["word/document.xml"] = ET.tostring(
            document, encoding="utf-8", xml_declaration=True
        )
        reviewed = tmp / "reviewed.docx"
        with ZipFile(reviewed, "w") as archive:
            for filename, data in entries.items():
                archive.writestr(filename, data)
        assert original_math == [
            ET.tostring(m) for m in document.findall(".//m:oMath", NS)
        ], "Word math changed"
        comparison = SimpleNamespace(
            reviewed=reviewed,
            baseline=baseline_dir / "manifest.json",
            output=tmp / "review.json",
        )
        saved_root = review.ROOT
        review.ROOT = tmp
        try:
            review.compare(comparison)
            report = json.loads(comparison.output.read_text())
            assert report["conflicts"] == 0
            assert len(report["comments"]) == 1
            assert intro_path in report["comments"][0]["source_candidates"]
            assert any(p["revisions"] for c in report["changes"] for p in c["after"])
            assert any(
                sentence in p["text"] for c in report["changes"] for p in c["after"]
            )
            assert any(
                "review verification" in p["text"]
                for c in report["changes"]
                for p in c["after"]
            )
            # Perform the agent's surgical merge into the original section snapshot.
            edited_intro = intro.replace("seventh", "7th", 1).replace(
                old_phrase, old_phrase + " review verification", 1
            )
            blocks = edited_intro.split("\n\n")
            anchor = next(i for i, block in enumerate(blocks) if old_phrase in block)
            blocks.insert(anchor + 1, sentence)
            edited_intro = (
                "% MM review 70003 [MM Test]: Please review this citation context.\n"
                + "\n\n".join(blocks)
            )
            (tmp / intro_path).write_text(edited_intro)
            assert re.findall(r"\\cite\w*(?:\[[^]]*\])*\{[^}]+\}", intro) == re.findall(
                r"\\cite\w*(?:\[[^]]*\])*\{[^}]+\}", edited_intro
            )
            review.compare(comparison)
            report = json.loads(comparison.output.read_text())
            assert report["conflicts"] > 0, "Later LaTeX edits must be reported"
        finally:
            review.ROOT = saved_root
        if args.compile:
            latex = tmp / "latex"
            shutil.copytree(ROOT / "build" / name / "latex", latex)
            main_tex = latex / "main.tex"
            original_main = main_tex.read_text()
            edited_main = original_main.replace("seventh", "7th", 1)
            assert old_phrase in edited_main
            edited_main = edited_main.replace(
                old_phrase, old_phrase + " review verification", 1
            )
            paragraph_end = edited_main.index("\n\n", edited_main.index(old_phrase))
            edited_main = (
                edited_main[:paragraph_end]
                + "\n\n"
                + sentence
                + edited_main[paragraph_end:]
            )
            assert re.findall(r"\\label\{[^}]+\}", original_main) == re.findall(
                r"\\label\{[^}]+\}", edited_main
            )
            assert re.findall(
                r"\\cite\w*(?:\[[^]]*\])*\{[^}]+\}", original_main
            ) == re.findall(r"\\cite\w*(?:\[[^]]*\])*\{[^}]+\}", edited_main)
            main_tex.write_text(edited_main)
            logfile = ROOT / "build" / f"mm-{args.manuscript}-return-test.log"
            with logfile.open("w") as log:
                for command in (
                    [
                        "pdflatex",
                        "-interaction=nonstopmode",
                        "-halt-on-error",
                        "main.tex",
                    ],
                    ["bibtex", "main"],
                    *[
                        [
                            "pdflatex",
                            "-interaction=nonstopmode",
                            "-halt-on-error",
                            "main.tex",
                        ]
                    ]
                    * 3,
                ):
                    subprocess.run(
                        command, cwd=latex, stdout=log, stderr=log, check=True
                    )
                checks = f"""
from pathlib import Path
from local.value_checks import extract_value_check_specs, run_value_checks
from local.{args.manuscript}_ghg_forcing_for_cmip7.value_checks import get_value_checks
source_dir = Path({str(tmp / "manuscripts" / name)!r})
generated_dirs = [Path({str(ROOT / "tables" / name)!r}),
                  Path({str(ROOT / "build" / name)!r})]
paths = list(source_dir.glob('*.tex'))
for directory in generated_dirs:
    paths.extend(directory.glob('*.tex'))
specs = [s for p in paths for s in extract_value_check_specs(p.read_text(), file=p)]
r = run_value_checks(specs, get_value_checks())
print(r.to_str())
r.raise_if_not_ok()
"""
                subprocess.run(
                    [str(ROOT / ".venv/bin/python"), "-c", checks],
                    cwd=ROOT,
                    stdout=log,
                    stderr=log,
                    check=True,
                )
            shutil.copy2(
                latex / "main.pdf",
                ROOT / "build" / f"mm-{args.manuscript}-return-test.pdf",
            )
    print(
        json.dumps(
            {
                "manuscript": name,
                "status": "passed",
                "compiled": args.compile,
                "checks": [
                    "tracked edit",
                    "ordinary edit",
                    "inserted paragraph",
                    "anchored comment",
                    "source conflict",
                    "citation preservation",
                    "native math preservation",
                ],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
