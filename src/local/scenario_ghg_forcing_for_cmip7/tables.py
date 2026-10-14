"""
Generation of the scenario GHG manuscript's tables
"""

from __future__ import annotations

from pathlib import Path

from loguru import logger

from local.scenario_ghg_forcing_for_cmip7.bundle import (
    ONE_BOX_RERUN_UNIT,
    load_inverse_emissions,
    load_one_box_rerun,
)

SUMMARY_TABLE_YEARS = (2023, 2050, 2100, 2200, 2300, 2500)
"""Years shown in the summary table"""

SUMMARY_TABLE_EMISSIONS_UNIT = "kt CCl4 / yr"
"""Unit in which the summary table shows the inverse emissions"""

SUMMARY_TABLE_COLUMNS = (
    ("concentration", "Concentration", ONE_BOX_RERUN_UNIT, "{:.1f}"),
    ("emissions", "Inverse emissions", SUMMARY_TABLE_EMISSIONS_UNIT, "{:.2f}"),
)
"""Columns of the summary table

Key, heading, unit and format of each column.
"""


def generate_summary_table(
    outfile: Path, inverse_emissions_file: Path, intermediate_file: Path
) -> Path:
    """
    Generate the summary table

    Only the tabular is generated,
    the caption etc. are in the manuscript's text.

    Parameters
    ----------
    outfile
        File in which to write the table

    inverse_emissions_file
        Inverse emissions written by the re-run inverse emissions notebook

    intermediate_file
        Intermediates saved by the re-run inverse emissions notebook

    Returns
    -------
    :
        `outfile`
    """
    years = list(SUMMARY_TABLE_YEARS)
    values = {
        "concentration": load_one_box_rerun(intermediate_file).loc[years].to_numpy(),
        "emissions": load_inverse_emissions(
            inverse_emissions_file, unit=SUMMARY_TABLE_EMISSIONS_UNIT
        )
        .loc[years]
        .to_numpy(),
    }

    rows = []
    for i, year in enumerate(SUMMARY_TABLE_YEARS):
        cells = [
            str(year),
            *(fmt.format(values[key][i]) for key, _, _, fmt in SUMMARY_TABLE_COLUMNS),
        ]
        rows.append(f"    {' & '.join(cells)} \\\\")

    header = " & ".join(
        ["Year", *(heading for _, heading, _, _ in SUMMARY_TABLE_COLUMNS)]
    )
    # Units written the way the text writes them,
    # so the manuscript's replacements format them
    sub_header = " & ".join(
        ["", *(f"[{unit}]" for _, _, unit, _ in SUMMARY_TABLE_COLUMNS)]
    ).replace("CCl4", "CCl_4")

    res = "\n".join(
        [
            r"\begin{tabular}{lrr}",
            r"\tophline",
            f"{header} \\\\",
            f"{sub_header} \\\\",
            r"\middlehline",
            *rows,
            r"\bottomhline",
            r"\end{tabular}",
            "",
        ]
    )

    outfile.parent.mkdir(exist_ok=True, parents=True)
    logger.info(f"Writing {outfile}")
    outfile.write_text(res)

    return outfile
