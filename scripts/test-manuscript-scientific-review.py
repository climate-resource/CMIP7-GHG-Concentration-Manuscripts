#!/usr/bin/env python3
"""Verify native revision semantics against real exported manuscript objects."""

import argparse
import copy
import importlib.util
import json
from pathlib import Path
from zipfile import ZipFile

from lxml import etree as E

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "scientific_review", ROOT / "scripts/manuscript-scientific-review.py"
)
review = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(review)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    args = parser.parse_args()
    with ZipFile(args.baseline.parent / "baseline.docx") as archive:
        parts = {name: archive.read(name) for name in archive.namelist()}
    reviewer = review.Reviewer(parts)
    # Add a disposable mixed-run paragraph with an equation between prose spans.
    paragraph = review.node("p")
    paragraph.append(reviewer.run(None, "The old "))
    paragraph.append(reviewer.run(None, "result is "))
    math = E.SubElement(paragraph, f"{{{review.M}}}oMath")
    math_run = E.SubElement(math, f"{{{review.M}}}r")
    E.SubElement(math_run, f"{{{review.M}}}t").text = "x"
    paragraph.append(reviewer.run(None, " beside (Smith, 2020)."))
    paragraph.append(reviewer.run(None, " See "))
    internal = review.node("hyperlink", anchor="fig:test")
    internal.append(reviewer.run(None, "70"))
    paragraph.append(internal)
    paragraph.append(reviewer.run(None, "b and "))
    external = review.node("hyperlink")
    external.set(
        "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id",
        "rIdTestOriginal",
    )
    external.append(reviewer.run(None, "https://doi.org/https://doi.org/example"))
    paragraph.append(external)
    # Insert the same original paragraph into the comparison package, ensuring
    # the protected-object check also covers this test equation.
    before = E.fromstring(parts["word/document.xml"])
    before.find("w:body", review.NS).insert(0, copy.deepcopy(paragraph))
    parts["word/document.xml"] = E.tostring(before)
    reviewer.root.find("w:body", review.NS).insert(0, paragraph)
    index = len(reviewer.paragraphs)
    reviewer.paragraphs.append(paragraph)
    reviewer.original.append(review.text_of(paragraph))
    reviewer.expected.append(review.text_of(paragraph))
    reviewer.replace(
        {
            "id": "cross-run",
            "paragraph_index": index,
            "old": "old result",
            "new": "revised estimate",
        }
    )
    reviewer.replace(
        {
            "id": "citation-adjacent",
            "paragraph_index": index,
            "old": "beside",
            "new": "consistent with",
        }
    )
    reviewer.insert(
        {
            "id": "new-paragraph",
            "paragraph_index": index,
            "new": "A proposed scientific caveat.",
        }
    )
    reviewer.comment(
        {
            "id": "comment",
            "paragraph_index": index,
            "anchor": "(Smith, 2020)",
            "text": "Check this primary source.",
        }
    )
    reviewer.replace(
        {
            "id": "doi-target",
            "paragraph_index": index,
            "old": "https://doi.org/https://doi.org/example",
            "new": "https://doi.org/example",
        }
    )
    reviewer.replace(
        {"id": "label-suffix", "paragraph_index": index, "old": "70b", "new": "C3b"}
    )
    reviewer.repair_cross_references({"fig:test": "C3"})
    assert review.text_of(internal) == "C3"
    state = E.tostring(paragraph)
    try:
        reviewer.replace(
            {
                "id": "unsafe",
                "paragraph_index": index,
                "old": "is x beside",
                "new": "wrong",
            }
        )
    except ValueError:
        pass
    else:
        raise AssertionError("Unsafe anchor was accepted")
    assert E.tostring(paragraph) == state
    # Protected math rejection is checked separately with a valid current anchor.
    try:
        reviewer.replace(
            {
                "id": "math",
                "paragraph_index": index,
                "old": "is x consistent",
                "new": "wrong",
            }
        )
    except ValueError as error:
        assert "protected equation" in str(error)
    else:
        raise AssertionError("Equation replacement was accepted")
    assert E.tostring(paragraph) == state
    result = reviewer.finish()
    assert result["comments"] == 1
    assert result["revision_elements"] >= 6
    assert review.text_of(paragraph) == (
        "The revised estimate is x consistent with (Smith, 2020). "
        "See C3b and https://doi.org/example"
    )
    assert not review.text_of(external)
    assert review.text_of(external, True) == "https://doi.org/https://doi.org/example"
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
