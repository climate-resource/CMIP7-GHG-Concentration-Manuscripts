"""
Generate GHG listing for use in markdown
"""

from pathlib import Path


def main():  # noqa: PLR0912, PLR0915
    """Generate the listing"""
    read_path = (
        Path(__file__).parents[2]
        / "CMIP-GHG-Concentration-Generation"
        / "output-bundles"
        / "v1.0.0"
        / "data"
        / "processed"
        / "esgf-ready"
        / "input4MIPs"
        / "CMIP7"
        / "CMIP"
        / "CR"
        / "CR-CMIP-1-0-0"
        / "atmos"
        / "mon"
    )

    big_three = []
    cfcs = []
    hcfcs = []
    halons = []
    other_ods = []
    hfcs = []
    pfcs = []
    other = []
    for d in sorted(read_path.iterdir()):
        out = d.name.upper()
        if "EQ" in out:
            continue

        for o, n in (
            ("CL", "Cl"),
            ("BR", "Br"),
            ("HALON", "Halon "),
            ("HFC", "HFC-"),
            ("CFC", "CFC-"),
            ("EA", "ea"),
            ("MEE", "mee"),
            ("FA", "fa"),
            ("134A", "134a"),
            ("143A", "143a"),
            ("152A", "152a"),
            ("142B", "142b"),
            ("141B", "141b"),
            ("MFC", "mfc"),
        ):
            out = out.replace(o, n)

        if not any(v in out for v in ("HFC", "HCFC", "Halon", "CFC")):
            for i in range(20)[::-1]:
                out = out.replace(str(i), f"<sub>{i}</sub>")
                out = (
                    out.replace("<sub><sub>", "<sub>")
                    .replace("</sub></sub>", "</sub>")
                    .replace("</sub><sub>", "")
                )

        if d.name in ("co2", "ch4", "n2o"):
            big_three.append(out)

        elif d.name.startswith("cfc"):
            cfcs.append(out)

        elif d.name.startswith("hcfc"):
            hcfcs.append(out)

        elif d.name.startswith("halon"):
            halons.append(out)

        elif d.name.startswith("ch") or d.name.startswith("ccl4"):
            other_ods.append(out)

        elif d.name.startswith("hfc"):
            hfcs.append(out)

        elif d.name.startswith("c"):
            pfcs.append(out)

        else:
            other.append(out)

    print(f"- major greenhouse gases ({len(big_three)})")
    print(f"    - {', '.join(big_three)}")
    print(
        f"- ozone-depleting substances ({len(cfcs) + len(hcfcs) + len(halons) + len(other_ods)})"  # noqa: E501
    )
    print(f"    - CFCs ({len(cfcs)})")
    print(f"        - {', '.join(cfcs)}")
    print(f"    - HCFCs ({len(hcfcs)})")
    print(f"        - {', '.join(hcfcs)}")
    print(f"    - Halons ({len(halons)})")
    print(f"        - {', '.join(halons)}")
    print(f"    - other ozone-depleting substances ({len(other_ods)})")
    print(f"        - {', '.join(other_ods)}")
    print(f"- ozone fluorinated compounds ({len(hfcs) + len(pfcs) + len(other)})")
    print(f"    - HFCs ({len(hfcs)})")
    print(f"        - {', '.join(hfcs)}")
    print(f"    - PFCs ({len(pfcs)})")
    print(f"        - {', '.join(pfcs)}")
    print(f"    - other ({len(other)})")
    print(f"        - {', '.join(other)}")


if __name__ == "__main__":
    main()
