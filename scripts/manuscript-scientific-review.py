#!/usr/bin/env python3
"""Apply checked scientific proposals to a saved Word export as native revisions.

Canonical LaTeX is never edited. Exact source proposals are emitted separately.
Run using the document artifact runtime (lxml required).
"""

from __future__ import annotations

import argparse
import copy
import difflib
import hashlib
import json
import re
import shutil
import uuid
from datetime import datetime, timezone
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

from lxml import etree as E

W = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
M = "http://schemas.openxmlformats.org/officeDocument/2006/math"
R = "http://schemas.openxmlformats.org/package/2006/relationships"
CT = "http://schemas.openxmlformats.org/package/2006/content-types"
NS = {"w": W, "m": M}
AUTHOR = "MM scientific review"


def tag(name):
    return f"{{{W}}}{name}"


def node(name, **attrs):
    el = E.Element(tag(name))
    for key, value in attrs.items():
        el.set(tag(key), str(value))
    return el


def text_of(el, original=False):
    if el.tag in (tag("ins"), tag("moveTo")) and original:
        return ""
    if el.tag in (tag("del"), tag("moveFrom")) and not original:
        return ""
    if el.tag in (tag("t"), tag("delText"), f"{{{M}}}t"):
        return el.text or ""
    if el.tag in (tag("tab"), tag("br")):
        return " "
    return "".join(text_of(child, original) for child in el)


def segments(el):
    """Accepted character positions, including protected math and breaks."""
    result = []
    offset = 0

    def visit(current):
        nonlocal offset
        if current.tag in (tag("del"), tag("moveFrom")):
            return
        if current.tag in (tag("t"), f"{{{M}}}t", tag("tab"), tag("br")):
            value = current.text or "" if current.tag.endswith("}t") else " "
            result.append((offset, offset + len(value), current))
            offset += len(value)
        else:
            for child in current:
                visit(child)

    visit(el)
    return result


class Reviewer:
    def __init__(self, parts):
        self.parts = parts
        self.root = E.fromstring(parts["word/document.xml"])
        self.paragraphs = self.root.findall(".//w:body//w:p", NS)
        self.original = [text_of(p) for p in self.paragraphs]
        self.counter = max(
            [0]
            + [int(v) for v in self.root.xpath("//@w:id", namespaces=NS) if v.isdigit()]
        )
        self.date = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        self.comments = (
            E.fromstring(parts["word/comments.xml"])
            if "word/comments.xml" in parts
            else node("comments")
        )
        self.comment_counter = max(
            [-1] + [int(c.get(tag("id"))) for c in self.comments]
        )
        self.applied = []
        self.expected = list(self.original)
        self.inserted = []

    def revision(self, kind):
        self.counter += 1
        return node(kind, id=self.counter, author=AUTHOR, date=self.date)

    def run(self, source_run, value, deleted=False):
        run = node("r")
        props = source_run.find("w:rPr", NS) if source_run is not None else None
        if props is not None:
            run.append(copy.deepcopy(props))
        content = node("delText" if deleted else "t")
        content.set("{http://www.w3.org/XML/1998/namespace}space", "preserve")
        content.text = value
        run.append(content)
        return run

    def add_text(self, parent, source_run, value):
        """Keep newly inserted chemical formulae and common units readable."""
        pattern = re.compile(
            r"\b(?:CO2|CH4|N2O|C8F18|CF4|C2F6|C3F8|NF3|SF6|SO2|O3|N2|O2|H2O|W/m2)\b"
        )
        position = 0
        for match in pattern.finditer(value):
            if match.start() > position:
                parent.append(self.run(source_run, value[position : match.start()]))
            for part in re.findall(r"\d+|\D+", match[0]):
                run = self.run(source_run, part)
                if part.isdigit():
                    props = run.find("w:rPr", NS)
                    if props is None:
                        props = node("rPr")
                        run.insert(0, props)
                    existing = props.find("w:vertAlign", NS)
                    if existing is not None:
                        props.remove(existing)
                    props.append(
                        node(
                            "vertAlign",
                            val="superscript" if match[0] == "W/m2" else "subscript",
                        )
                    )
                parent.append(run)
            position = match.end()
        if position < len(value):
            parent.append(self.run(source_run, value[position:]))

    def edit_segment(self, text_node, start, end, replacement):
        """Split one text run, preserving all unaffected nodes and run properties."""
        run = text_node.getparent()
        if run.tag != tag("r"):
            raise ValueError("Text is not in a normal Word run")
        if any(a.tag in (tag("ins"), tag("del")) for a in run.iterancestors()):
            raise ValueError("Overlapping proposal or pre-existing revision")
        parent = run.getparent()
        if parent.tag == tag("hyperlink") and not parent.get(tag("anchor")):
            raise ValueError("Hyperlink text requires individual relationship review")
        value = text_node.text or ""
        children = list(run)
        position = children.index(text_node)
        before = node("r")
        after = node("r")
        props = run.find("w:rPr", NS)
        for target in (before, after):
            if props is not None:
                target.append(copy.deepcopy(props))
        # Move original nodes so previously indexed text-node objects stay valid.
        for child in children[:position]:
            if child.tag != tag("rPr"):
                before.append(child)
        if start:
            before.append(self.run(None, value[:start])[0])
        if end < len(value):
            text_node.text = value[end:]
            after.append(text_node)
        for child in children[position + 1 :]:
            after.append(child)
        items = []
        if len(before) > (1 if props is not None else 0):
            items.append(before)
        if end > start:
            deletion = self.revision("del")
            deletion.append(self.run(run, value[start:end], True))
            items.append(deletion)
        if replacement:
            insertion = self.revision("ins")
            self.add_text(insertion, run, replacement)
            items.append(insertion)
        if len(after) > (1 if props is not None else 0):
            items.append(after)
        index = parent.index(run)
        parent.remove(run)
        for item in items:
            parent.insert(index, item)
            index += 1

    def change(self, paragraph, start, end, replacement):
        spans = segments(paragraph)
        if start == end:
            matches = [
                (a, b, n) for a, b, n in spans if n.tag == tag("t") and a <= start <= b
            ]
            if not matches:
                raise ValueError("No safe text run at insertion boundary")
            a, _, n = matches[-1]
            self.edit_segment(n, start - a, start - a, replacement)
            return
        matches = [(a, b, n) for a, b, n in spans if b > start and a < end]
        if not matches or any(n.tag != tag("t") for _, _, n in matches):
            raise ValueError(
                "Proposed replacement intersects protected equation or break"
            )
        if any(
            n.getparent().getparent().tag == tag("hyperlink")
            and not n.getparent().getparent().get(tag("anchor"))
            for _, _, n in matches
        ):
            raise ValueError("Proposed replacement intersects a hyperlink")
        for index in range(len(matches) - 1, -1, -1):
            a, b, n = matches[index]
            self.edit_segment(
                n, max(0, start - a), min(b, end) - a, replacement if index == 0 else ""
            )

    def replace(self, proposal):
        index = proposal["paragraph_index"]
        p = self.paragraphs[index]
        old, new = proposal["old"], proposal["new"]
        current = text_of(p)
        if not old or current.count(old) != 1:
            raise ValueError(
                f"{proposal['id']}: old text not unique in paragraph {index}"
            )
        offset = current.index(old)
        # A corrected DOI must not inherit the old erroneous hyperlink target.
        # Retain the old linked text as a deletion and insert the corrected URL
        # outside that link. Rejecting restores both original text and target.
        links = [
            link
            for link in p.findall(".//w:hyperlink", NS)
            if not link.get(tag("anchor")) and text_of(link) == old
        ]
        if len(links) == 1 and re.fullmatch(r"https://\S+", new):
            link = links[0]
            if any(c.tag != tag("r") for c in link):
                raise ValueError("Complex external hyperlink needs individual review")
            deletion = self.revision("del")
            for child in list(link):
                for t in child.iter(tag("t")):
                    t.tag = tag("delText")
                deletion.append(child)
            link.append(deletion)
            insertion = self.revision("ins")
            self.add_text(insertion, None, new)
            link.addnext(insertion)
            self.expected[index] = current[:offset] + new + current[offset + len(old) :]
            self.applied.append(proposal["id"])
            return
        old_tokens = re.findall(r"\w+|\s+|[^\w\s]", old)
        new_tokens = re.findall(r"\w+|\s+|[^\w\s]", new)
        boundaries = [0]
        for token in old_tokens:
            boundaries.append(boundaries[-1] + len(token))
        matcher = difflib.SequenceMatcher(None, old_tokens, new_tokens, autojunk=False)
        for kind, a, b, c, d in reversed(matcher.get_opcodes()):
            if kind != "equal":
                self.change(
                    p,
                    offset + boundaries[a],
                    offset + boundaries[b],
                    "".join(new_tokens[c:d]),
                )
        expected = current[:offset] + new + current[offset + len(old) :]
        if text_of(p) != expected:
            raise AssertionError(f"Accepted text mismatch: {proposal['id']}")
        self.expected[index] = expected
        self.applied.append(proposal["id"])

    def insert(self, proposal, anchor=None, bibliography=False):
        anchor = (
            self.paragraphs[proposal["paragraph_index"]] if anchor is None else anchor
        )
        paragraph = node("p")
        props = node("pPr")
        props.append(node("pStyle", val="Bibliography" if bibliography else "BodyText"))
        mark = node("rPr")
        mark.append(self.revision("ins"))
        props.append(mark)
        paragraph.append(props)
        insertion = self.revision("ins")
        self.add_text(insertion, None, proposal["new"])
        paragraph.append(insertion)
        anchor.addnext(paragraph)
        self.inserted.append(paragraph)
        self.applied.append(proposal["id"])
        return paragraph

    def comment(self, proposal):
        paragraph = self.paragraphs[proposal["paragraph_index"]]
        self.comment_counter += 1
        cid = self.comment_counter
        comment = node("comment", id=cid, author=AUTHOR, date=self.date, initials="MMR")
        para = node("p")
        quote = proposal.get("anchor", "")
        value = proposal["text"]
        if quote:
            if (
                quote not in text_of(paragraph)
                and quote not in self.original[proposal["paragraph_index"]]
            ):
                raise ValueError(f"{proposal['id']}: comment anchor not found")
            value = f"On “{quote}”: {value}"
        evidence = proposal.get("evidence", [])
        urls = list(dict.fromkeys(e["url"] for e in evidence if e.get("url")))
        if urls:
            value += " Sources: " + "; ".join(urls)
        para.append(self.run(None, value))
        comment.append(para)
        self.comments.append(comment)
        start = node("commentRangeStart", id=cid)
        paragraph.insert(1 if paragraph.find("w:pPr", NS) is not None else 0, start)
        paragraph.append(node("commentRangeEnd", id=cid))
        reference = node("r")
        reference.append(node("commentReference", id=cid))
        paragraph.append(reference)

    def finish(self):
        for index, paragraph in enumerate(self.paragraphs):
            if text_of(paragraph, True) != self.original[index]:
                raise AssertionError(
                    f"Rejecting revisions does not restore paragraph {index}"
                )
            if text_of(paragraph) != self.expected[index]:
                raise AssertionError(f"Unexpected accepted text in paragraph {index}")
        if any(text_of(p, True) for p in self.inserted):
            raise AssertionError("Inserted paragraph survives rejected text")
        before = E.fromstring(self.parts["word/document.xml"])
        for xpath in (".//m:oMath", ".//w:drawing", ".//w:bookmarkStart"):
            if [E.tostring(n) for n in before.findall(xpath, NS)] != [
                E.tostring(n) for n in self.root.findall(xpath, NS)
            ]:
                raise AssertionError(f"Protected objects changed: {xpath}")
        revisions = self.root.xpath("//w:ins|//w:del", namespaces=NS)
        ids = [r.get(tag("id")) for r in revisions]
        if len(ids) != len(set(ids)):
            raise AssertionError("Duplicate revision identifiers")
        self.parts["word/document.xml"] = E.tostring(
            self.root, xml_declaration=True, encoding="UTF-8", standalone=True
        )
        settings = E.fromstring(self.parts["word/settings.xml"])
        if settings.find("w:trackRevisions", NS) is None:
            settings.append(node("trackRevisions"))
        self.parts["word/settings.xml"] = E.tostring(
            settings, xml_declaration=True, encoding="UTF-8", standalone=True
        )
        if len(self.comments):
            self.parts["word/comments.xml"] = E.tostring(
                self.comments, xml_declaration=True, encoding="UTF-8", standalone=True
            )
            rels = E.fromstring(self.parts["word/_rels/document.xml.rels"])
            if not any(r.get("Type", "").endswith("/comments") for r in rels):
                E.SubElement(
                    rels,
                    f"{{{R}}}Relationship",
                    Id="rIdMMScientificComments",
                    Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/comments",
                    Target="comments.xml",
                )
            self.parts["word/_rels/document.xml.rels"] = E.tostring(
                rels, xml_declaration=True, encoding="UTF-8", standalone=True
            )
            types = E.fromstring(self.parts["[Content_Types].xml"])
            if not any(t.get("PartName") == "/word/comments.xml" for t in types):
                E.SubElement(
                    types,
                    f"{{{CT}}}Override",
                    PartName="/word/comments.xml",
                    ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.comments+xml",
                )
            self.parts["[Content_Types].xml"] = E.tostring(
                types, xml_declaration=True, encoding="UTF-8", standalone=True
            )
        return {
            "applied_edits": self.applied,
            "revision_elements": len(revisions),
            "comments": len(self.comments),
            "original_text_restored": True,
            "equations_figures_bookmarks_preserved": True,
        }

    def repair_cross_references(self, labels):
        """Track export-only fixes against authoritative compiled LaTeX labels."""
        for index, paragraph in enumerate(self.paragraphs):
            for link in paragraph.findall(".//w:hyperlink", NS):
                anchor = link.get(tag("anchor"))
                expected = labels.get(anchor)
                current = text_of(link)
                if expected is None or expected == current:
                    continue
                span = [
                    (a, b)
                    for a, b, n in segments(paragraph)
                    if link in n.iterancestors()
                ]
                start, end = span[0][0], span[-1][1]
                original = text_of(paragraph)
                self.change(paragraph, start, end, expected)
                self.expected[index] = original[:start] + expected + original[end:]
                self.applied.append(f"export-reference:{index}:{anchor}")


def proposed_latex(proposals, baseline, folder):
    texts = {p: item["text"] for p, item in baseline["snapshots"].items()}
    statuses = []
    for edit in proposals.get("edits", []):
        path, old, new = (
            edit.get(k) for k in ("source_file", "latex_before", "latex_after")
        )
        status = "requires individual LaTeX transfer"
        if path in texts and old and new is not None:
            count = texts[path].count(old)
            if count == 1:
                texts[path] = texts[path].replace(old, new, 1)
                status = "proposed source copy updated"
            else:
                status = f"source anchor occurs {count} times; manual review required"
        statuses.append({"id": edit["id"], "status": status})
    bibliography = "references/references.bib"
    for reference in proposals.get("new_references", []):
        key = reference["key"]
        if re.search(r"@\w+\s*\{\s*" + re.escape(key) + r"\s*,", texts[bibliography]):
            raise ValueError(f"New reference key already exists: {key}")
        texts[bibliography] += "\n\n" + reference["bibtex"] + "\n"
    patch = []
    for path, text in texts.items():
        old = baseline["snapshots"][path]["text"]
        if text != old:
            target = folder / "proposed-latex" / path
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_text(text)
            patch.extend(
                difflib.unified_diff(
                    old.splitlines(True),
                    text.splitlines(True),
                    fromfile="a/" + path,
                    tofile="b/" + path,
                )
            )
    (folder / "manuscript_track.patch").write_text("".join(patch))
    (folder / "latex-transfer-status.json").write_text(json.dumps(statuses, indent=2))
    return statuses


def write_review_report(proposals, bundle, qa):
    lines = [
        "# Scientific review proposals",
        "",
        str(proposals.get("summary", "")),
        "",
        f"Native Word revision elements: {qa['revision_elements']}; "
        f"anchored comments: {qa['comments']}.",
        "",
        "All revisions are proposals. Numerical analyses requested in comments "
        "have not been performed. Canonical LaTeX is unchanged.",
        "",
    ]
    for group, title in (
        ("edits", "Tracked edits"),
        ("comments", "Open review comments"),
    ):
        lines += ["## " + title, ""]
        for item in proposals.get(group, []):
            lines += [
                f"### {item['id']} ({item.get('severity', 'review')})",
                "",
                f"Baseline paragraph {item['paragraph_index']}; "
                f"source: {item.get('source_file', 'see paragraph map')}.",
                "",
                item.get("rationale", item.get("text", "")),
                "",
            ]
            if group == "edits":
                lines += ["Proposed text: " + item["new"], ""]
            for evidence in item.get("evidence", []):
                lines += [
                    f"- [{evidence.get('title', 'Primary source')}]({evidence['url']})",
                    "",
                ]
    lines += [
        "## Reference audit coverage",
        "",
        json.dumps(proposals.get("coverage", {}), ensure_ascii=False, indent=2),
        "",
    ]
    (bundle / "REVIEW.md").write_text("\n".join(lines))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("proposals", type=Path)
    parser.add_argument("--baseline", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    proposals = json.loads(args.proposals.read_text())
    baseline = json.loads(args.baseline.read_text())
    if proposals["baseline_export_id"] != baseline["export_id"]:
        raise ValueError("Proposal and baseline export identifiers differ")
    if args.output.exists():
        raise FileExistsError("Refusing to overwrite an existing review draft")
    with ZipFile(args.baseline.parent / "baseline.docx") as archive:
        parts = {name: archive.read(name) for name in archive.namelist()}
    reviewer = Reviewer(parts)
    # Reverse index/order prevents consecutive paragraph insertions being inverted.
    for edit in reversed(proposals.get("edits", [])):
        if edit["kind"] == "replace":
            reviewer.replace(edit)
        elif edit["kind"] == "insert_after":
            reviewer.insert(edit)
        else:
            raise ValueError(f"Unknown edit kind: {edit['kind']}")
    reviewer.repair_cross_references(baseline["labels"])
    for comment in proposals.get("comments", []):
        reviewer.comment(comment)
    # Suggested references are genuine tracked additions, anchored after the
    # existing bibliography; this avoids silently renumbering existing citations.
    refs = proposals.get("new_references", [])
    if refs:
        heading = next(
            i
            for i, text in enumerate(reviewer.original)
            if text.strip() == "References"
        )
        candidates = [
            p
            for p in baseline["paragraphs"]
            if p["index"] > heading and p.get("source") == "references/references.bib"
        ]
        # Original export ends the references at the review-only source TODOs.
        todo = next(
            (
                i
                for i, text in enumerate(reviewer.original)
                if i > heading and "Source manuscript TODO" in text
            ),
            len(reviewer.original),
        )
        last = max(p["index"] for p in candidates if p["index"] < todo)
        anchor = reviewer.paragraphs[last]
        for ref in refs:
            anchor = reviewer.insert(
                {
                    "id": "reference:" + ref["key"],
                    "paragraph_index": last,
                    "new": ref["word_reference"],
                },
                anchor,
                True,
            )
    qa = reviewer.finish()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with ZipFile(args.output, "w", ZIP_DEFLATED) as archive:
        for name, data in reviewer.parts.items():
            archive.writestr(name, data)
    review_id = (
        datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        + "-"
        + uuid.uuid4().hex[:8]
    )
    bundle = args.baseline.parent.parent / (
        baseline["export_id"] + "-trackdraft-" + review_id
    )
    bundle.mkdir()
    shutil.copy2(args.output, bundle / "baseline.docx")
    shutil.copy2(args.proposals, bundle / "proposals.json")
    audit = args.proposals.parent / "reference-audit.json"
    if audit.exists():
        shutil.copy2(audit, bundle / audit.name)
    statuses = proposed_latex(proposals, baseline, bundle)
    for ref in refs:
        with (bundle / "suggested-references.bib").open("a") as handle:
            handle.write(ref["bibtex"] + "\n\n")
    report = {
        "review_id": review_id,
        "parent_export_id": baseline["export_id"],
        "source_commit": baseline["source_commit"],
        "source_baseline": str(args.baseline.resolve()),
        "output": str(args.output.resolve()),
        "sha256": hashlib.sha256(args.output.read_bytes()).hexdigest(),
        "qa": qa,
        "latex_transfer": statuses,
        "note": (
            "Unaccepted scientific suggestions. Canonical LaTeX unchanged. "
            "Resolve or explicitly select revisions before assisted transfer."
        ),
    }
    (bundle / "review-manifest.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2)
    )
    write_review_report(proposals, bundle, qa)
    print(
        json.dumps(
            {"output": str(args.output), "review_bundle": str(bundle), "qa": qa},
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
