"""
Generate GHG listing for use in markdown

Copied from `generate-ghg-listing.py` for the historical EO-update guide
(CO2 and CH4, extended with satellite data), which reads from a separate
output bundle to the original historical dataset.
"""

from pathlib import Path


def main():
    """Generate the listing"""
    # TODO: point this at the final, published EO-updated dataset once it
    # exists. For now, this points at a local output bundle from
    # CMIP-GHG-Concentration-Generation (eo-update, CR-CMIP7Plus).
    read_path = (
        Path(__file__).parents[2]
        / "CMIP-GHG-Concentration-Generation"
        / "output-bundles"
        / "eo-update"
        / "data"
        / "processed"
        / "esgf-ready"
        / "EO-update-input4MIPs"
        / "CMIP7Plus"
        / "CMIP"
        / "CR"
        / "CR-CMIP7Plus"
        / "atmos"
        / "mon"
    )

    species = []
    for d in sorted(read_path.iterdir()):
        out = d.name.upper()
        for i in range(10)[::-1]:
            out = out.replace(str(i), f"{{raw-latex}}`\\textsubscript{{{i}}}`")

        species.append(out)

    print(f"- major greenhouse gases ({len(species)})")
    print(f"    - {', '.join(species)}")


if __name__ == "__main__":
    main()
