"""
Compile the tex inputs for the historical GHG manuscript
"""

import json
from pathlib import Path
from typing import Annotated

import typer

from local.cmip_ghg_generation import (
    DEFAULT_BUNDLE_DIR,
    DEFAULT_ORIGINAL_RUN_NOTEBOOKS_DIR,
)
from local.historical_ghg_forcing_for_cmip7 import (
    C4F10_LIKE_GASES,
    CFC12_LIKE_GASES,
    generate_c4f10_like_methods_figure,
    generate_c8f18_methods_figure,
    generate_cfc12_like_methods_figure,
    generate_cfc12_like_obs_network_sources_list,
    generate_cfc12_like_per_gas_table,
    generate_ch4_methods_figure,
    generate_co2_methods_figure,
    generate_n2o_methods_figure,
    generate_results_figure_for_gas,
)
from local.historical_ghg_forcing_for_cmip7.comparison_data import (
    RADIATIVE_EFFICIENCIES,
)


def parse_figure_file(
    value: str, gases: tuple[str, ...], group: str
) -> tuple[str, Path]:
    """
    Parse a gas and a file to write that gas' figure in

    Parameters
    ----------
    value
        Value to parse, as `"<gas>=<file>"`

    gases
        Gases which the figure this value is for can be drawn for

    group
        How the group of gases is named in an error message

    Returns
    -------
    :
        The gas and the file to write its figure in

    Raises
    ------
    typer.BadParameter
        `value` is not a gas and a file,
        or the gas is not one of `group`
    """
    gas, _, figure_file = value.partition("=")
    if not figure_file:
        msg = f"Expected '<gas>=<file>', got {value!r}"
        raise typer.BadParameter(msg)

    if gas not in gases:
        msg = f"{gas!r} is not one of {group}. Expected one of: {', '.join(gases)}"
        raise typer.BadParameter(msg)

    return gas, Path(figure_file)


def main(  # noqa: PLR0913
    co2_methods_figure_file: Annotated[
        Path,
        typer.Option(
            help="Path to in which to write the CO2 methods figure. ",
            dir_okay=False,
            file_okay=True,
        ),
    ],
    co2_methods_appendix_figure_file: Annotated[
        Path,
        typer.Option(
            help=(
                "Path to in which to write the CO2 methods appendix figure, "
                "which shows the working behind the methods figure."
            ),
            dir_okay=False,
            file_okay=True,
        ),
    ],
    ch4_methods_figure_file: Annotated[
        Path,
        typer.Option(
            help="Path to in which to write the CH4 methods figure. ",
            dir_okay=False,
            file_okay=True,
        ),
    ],
    ch4_methods_appendix_figure_file: Annotated[
        Path,
        typer.Option(
            help=(
                "Path to in which to write the CH4 methods appendix figure, "
                "which shows the working behind the methods figure."
            ),
            dir_okay=False,
            file_okay=True,
        ),
    ],
    n2o_methods_figure_file: Annotated[
        Path,
        typer.Option(
            help="Path to in which to write the N2O methods figure. ",
            dir_okay=False,
            file_okay=True,
        ),
    ],
    n2o_methods_appendix_figure_file: Annotated[
        Path,
        typer.Option(
            help=(
                "Path to in which to write the N2O methods appendix figure, "
                "which shows the working behind the methods figure."
            ),
            dir_okay=False,
            file_okay=True,
        ),
    ],
    c8f18_methods_figure_file: Annotated[
        Path,
        typer.Option(
            help="Path to in which to write the C8F18 methods figure. ",
            dir_okay=False,
            file_okay=True,
        ),
    ],
    cfc12_like_obs_network_sources_list_file: Annotated[
        Path,
        typer.Option(
            help=(
                "Path in which to write the list of which observation networks "
                "observe which of the gases processed like CFC-12."
            ),
            dir_okay=False,
            file_okay=True,
        ),
    ],
    cfc12_like_per_gas_table_file: Annotated[
        Path,
        typer.Option(
            help=(
                "Path in which to write the table of per-gas inputs and choices "
                "for the gases processed like CFC-12."
            ),
            dir_okay=False,
            file_okay=True,
        ),
    ],
    tex_inputs_manifest_file: Annotated[
        Path,
        typer.Option(
            help=(
                "Path in which to write the manifest of what was generated. "
                "It maps the tag each figure and table replaces in the latex "
                "to the file it was written in, "
                "and is passed straight on to the compilation script."
            ),
            dir_okay=False,
            file_okay=True,
        ),
    ],
    cfc12_like_methods_figure_file: Annotated[
        list[str] | None,
        typer.Option(
            help=(
                "Gas which is processed like CFC-12 "
                "and the path in which to write its methods figure, "
                "as '<gas>=<file>' e.g. 'cfc12=figures/cfc12_methods.pdf'. "
                "Repeat the option for each gas you want a figure for."
            ),
        ),
    ] = None,
    cfc12_like_methods_appendix_figure_file: Annotated[
        list[str] | None,
        typer.Option(
            help=(
                "Gas which is processed like CFC-12 "
                "and the path in which to write its methods appendix figure, "
                "as '<gas>=<file>'. "
                "Every gas given a methods figure must be given one of these too."
            ),
        ),
    ] = None,
    c4f10_like_methods_figure_file: Annotated[
        list[str] | None,
        typer.Option(
            help=(
                "Gas which is processed like C4F10 "
                "and the path in which to write its methods figure, "
                "as '<gas>=<file>' e.g. 'c4f10=figures/c4f10_methods.pdf'. "
                "Repeat the option for each gas you want a figure for."
            ),
        ),
    ] = None,
    results_figure_file: Annotated[
        list[str] | None,
        typer.Option(
            help=(
                "Gas and the path in which to write its results figure, "
                "as '<gas>=<file>' e.g. 'ch4=figures/ch4_results.pdf'. "
                "Repeat the option for each gas you want a figure for."
            ),
        ),
    ] = None,
    bundle_dir: Annotated[
        Path,
        typer.Option(
            help=("Directory in which to keep the original run's bundle."),
            dir_okay=True,
            file_okay=False,
        ),
    ] = DEFAULT_BUNDLE_DIR,
    original_run_notebooks_dir: Annotated[
        Path,
        typer.Option(
            help=(
                "Path to the original run's executed notebooks. "
                "Only needed if this repository doesn't already have a copy "
                "of the notebooks we re-run."
            ),
            dir_okay=True,
            file_okay=False,
        ),
    ] = DEFAULT_ORIGINAL_RUN_NOTEBOOKS_DIR,
    force_rerun: Annotated[
        bool,
        typer.Option(
            help=(
                "Re-run the original run's notebooks, "
                "even if we already have the data they produce."
            ),
        ),
    ] = False,
) -> None:
    """
    Create the inputs
    """
    # Parsed before anything is drawn, so a typo in a gas name
    # is caught now rather than three figures from now.
    cfc12_like_figures = [
        parse_figure_file(value, CFC12_LIKE_GASES, "the gases processed like CFC-12")
        for value in cfc12_like_methods_figure_file or []
    ]
    cfc12_like_appendix_figures = dict(
        parse_figure_file(value, CFC12_LIKE_GASES, "the gases processed like CFC-12")
        for value in cfc12_like_methods_appendix_figure_file or []
    )
    # The two figures are drawn together, so a gas needs both or neither
    if set(cfc12_like_appendix_figures) != {gas for gas, _ in cfc12_like_figures}:
        msg = (
            "Every gas processed like CFC-12 needs both a methods figure file "
            "and a methods appendix figure file. "
            f"Methods: {sorted(gas for gas, _ in cfc12_like_figures)}. "
            f"Appendix: {sorted(cfc12_like_appendix_figures)}."
        )
        raise typer.BadParameter(msg)

    c4f10_like_figures = [
        parse_figure_file(value, C4F10_LIKE_GASES, "the gases processed like C4F10")
        for value in c4f10_like_methods_figure_file or []
    ]
    # Every gas (and equivalent species) has a radiative efficiency,
    # so this is the list of gases.
    results_figures = [
        parse_figure_file(
            value,
            tuple(RADIATIVE_EFFICIENCIES),
            "the gases we produce",
        )
        for value in results_figure_file or []
    ]

    # The tag each output replaces in the latex, and the file it is written in.
    # The tags follow one convention, set here,
    # so the build script only has to say where to write each file.
    figures: dict[str, Path] = {}
    tables: dict[str, Path] = {}

    generate_co2_methods_figure(
        co2_methods_figure_file,
        co2_methods_appendix_figure_file,
        bundle_dir=bundle_dir,
        original_run_notebooks_dir=original_run_notebooks_dir,
        force_rerun=force_rerun,
    )
    figures["<co2-methods-figure>"] = co2_methods_figure_file
    figures["<co2-methods-appendix-figure>"] = co2_methods_appendix_figure_file

    generate_ch4_methods_figure(
        ch4_methods_figure_file,
        ch4_methods_appendix_figure_file,
        bundle_dir=bundle_dir,
        original_run_notebooks_dir=original_run_notebooks_dir,
        force_rerun=force_rerun,
    )
    figures["<ch4-methods-figure>"] = ch4_methods_figure_file
    figures["<ch4-methods-appendix-figure>"] = ch4_methods_appendix_figure_file

    generate_n2o_methods_figure(
        n2o_methods_figure_file,
        n2o_methods_appendix_figure_file,
        bundle_dir=bundle_dir,
        original_run_notebooks_dir=original_run_notebooks_dir,
        force_rerun=force_rerun,
    )
    figures["<n2o-methods-figure>"] = n2o_methods_figure_file
    figures["<n2o-methods-appendix-figure>"] = n2o_methods_appendix_figure_file

    # One figure per gas asked for: the gases processed like CFC-12
    # all share a figure, but there are thirty-four of them,
    # so which ones we draw is the caller's call.
    for gas, figure_file in cfc12_like_figures:
        generate_cfc12_like_methods_figure(
            gas,
            figure_file,
            cfc12_like_appendix_figures[gas],
            bundle_dir=bundle_dir,
            original_run_notebooks_dir=original_run_notebooks_dir,
            force_rerun=force_rerun,
        )
        figures[f"<{gas}-methods-figure>"] = figure_file
        figures[f"<{gas}-methods-appendix-figure>"] = cfc12_like_appendix_figures[gas]

    # The C4F10-like gases and C8F18 are built entirely from data
    # the original run left in the bundle,
    # so nothing here has a notebook to re-run.
    for gas, figure_file in c4f10_like_figures:
        generate_c4f10_like_methods_figure(
            gas,
            figure_file,
            bundle_dir=bundle_dir,
            force_rerun=force_rerun,
        )
        figures[f"<{gas}-methods-figure>"] = figure_file

    generate_c8f18_methods_figure(
        c8f18_methods_figure_file,
        bundle_dir=bundle_dir,
        force_rerun=force_rerun,
    )
    figures["<c8f18-methods-figure>"] = c8f18_methods_figure_file

    for gas, figure_file in results_figures:
        generate_results_figure_for_gas(
            gas,
            figure_file,
            bundle_dir=bundle_dir,
            original_run_notebooks_dir=original_run_notebooks_dir,
            force_rerun=force_rerun,
        )
        figures[f"<{gas}-results-figure>"] = figure_file

    # These summarise every gas processed like CFC-12,
    # so unlike the figures, they don't depend on which gases were asked for.
    generate_cfc12_like_obs_network_sources_list(
        cfc12_like_obs_network_sources_list_file,
        bundle_dir=bundle_dir,
        force_rerun=force_rerun,
    )
    tables["<cfc12-like-obs-network-sources-list>"] = (
        cfc12_like_obs_network_sources_list_file
    )

    generate_cfc12_like_per_gas_table(
        cfc12_like_per_gas_table_file,
        bundle_dir=bundle_dir,
        force_rerun=force_rerun,
    )
    tables["<cfc12-like-per-gas-table>"] = cfc12_like_per_gas_table_file

    # Absolute paths, so the manifest doesn't depend on
    # the directory the compilation script is run from.
    tex_inputs_manifest_file.parent.mkdir(parents=True, exist_ok=True)
    tex_inputs_manifest_file.write_text(
        json.dumps(
            {
                "figures": {tag: str(fp.resolve()) for tag, fp in figures.items()},
                "tables": {tag: str(fp.resolve()) for tag, fp in tables.items()},
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    typer.run(main)
