"""
Tasks which generate the figures and other tex inputs
"""

from pathlib import Path
from typing import Annotated

from pytask import Product, task

import local.historical_ghg_forcing_for_cmip7.author_contributions
import local.historical_ghg_forcing_for_cmip7.comparison_data
import local.historical_ghg_forcing_for_cmip7.equivalent_species
import local.scenario_ghg_forcing_for_cmip7.appendix_figures
import local.scenario_ghg_forcing_for_cmip7.emissions_figure
import local.scenario_ghg_forcing_for_cmip7.esgf_data
import local.scenario_ghg_forcing_for_cmip7.methods_data
import local.scenario_ghg_forcing_for_cmip7.methods_figures
import local.scenario_ghg_forcing_for_cmip7.results_figure
import local.scenario_ghg_forcing_for_cmip7.scenarios
from local.historical_ghg_forcing_for_cmip7.author_contributions import (
    generate_author_contributions,
)
from local.scenario_ghg_forcing_for_cmip7.appendix_figures import (
    APPENDIX_FIGURE_GROUPS,
    CFC12_LIKE_GASES,
    EMISSIONS_FIGURE_SPECIES,
    HFC125_LIKE_GASES,
    MAGICC_MISMATCH_GASES,
    N2O_LIKE_GASES,
    REGRESSION_GASES,
    RESULTS_FIGURE_GASES,
    generate_appendix_figures,
    get_tags_in_texts,
)
from local.scenario_ghg_forcing_for_cmip7.bundle import BUNDLE_DIR
from local.scenario_ghg_forcing_for_cmip7.config import (
    APPENDICES,
    APPENDIX_FIGURES_FILES,
    AUTHOR_CONTRIBUTIONS,
    BUNDLE_EXTRACTED_MARKERS,
    HISTORICAL_BUNDLE_EXTRACTED_MARKERS,
    MAGICC_SUMMARY,
    METADATA,
    METHODS_FIGURE_GASES,
    RERUN_CO2_SEASONALITY_FILES,
    RERUN_CO2_SEASONALITY_VL_CF_FILES,
    RERUN_LAT_GRADIENT_GASES,
    SECTIONS,
    WESTERN_ET_AL_2024_FIGURE,
    get_emissions_figure_file,
    get_magicc_mismatch_figure_file,
    get_methods_figure_file,
    get_regressions_figure_file,
    get_rerun_lat_gradient_files,
    get_results_figure_file,
)
from local.scenario_ghg_forcing_for_cmip7.emissions_figure import (
    CMIP6_EMISSIONS_FILE,
    generate_emissions_figure,
)
from local.scenario_ghg_forcing_for_cmip7.methods_figures import (
    generate_cfc12_like_figure,
    generate_ch4_figure,
    generate_co2_figure,
    generate_hfc125_like_figure,
    generate_magicc_mismatch_figure,
    generate_n2o_like_figure,
    generate_regressions_figure,
    generate_western_et_al_2024_figure,
)
from local.scenario_ghg_forcing_for_cmip7.results_figure import (
    generate_results_figure,
)

FIGURES_SOURCES = [
    Path(local.scenario_ghg_forcing_for_cmip7.methods_figures.__file__),
    Path(local.scenario_ghg_forcing_for_cmip7.methods_data.__file__),
]
"""The modules which make the figures (they are re-made if these change)"""

RESULTS_FIGURES_SOURCES = [
    Path(local.scenario_ghg_forcing_for_cmip7.results_figure.__file__),
    Path(local.scenario_ghg_forcing_for_cmip7.esgf_data.__file__),
    Path(local.scenario_ghg_forcing_for_cmip7.scenarios.__file__),
    # IGCC's record
    Path(local.historical_ghg_forcing_for_cmip7.comparison_data.__file__),
    Path(local.historical_ghg_forcing_for_cmip7.equivalent_species.__file__),
]
"""The modules which make the results figures (they are re-made if these change)"""

EMISSIONS_FIGURES_SOURCES = [
    *RESULTS_FIGURES_SOURCES,
    Path(local.scenario_ghg_forcing_for_cmip7.emissions_figure.__file__),
    Path(local.scenario_ghg_forcing_for_cmip7.methods_data.__file__),
]
"""The modules which make the emissions figures (they are re-made if these change)

They share the results figures' yearly panels.
"""


def get_rerun_files(gas: str) -> list[Path]:
    """Get the re-run outputs a gas' methods figure needs"""
    res = get_rerun_lat_gradient_files(gas) if gas in RERUN_LAT_GRADIENT_GASES else []
    if gas == "co2":
        res = [*res, *RERUN_CO2_SEASONALITY_FILES, *RERUN_CO2_SEASONALITY_VL_CF_FILES]

    return res


for _gas in METHODS_FIGURE_GASES:

    @task(id=_gas)
    def task_methods_figure(  # noqa: PLR0913
        gas: str = _gas,
        magicc_summary: Path = MAGICC_SUMMARY,
        extracted_markers: list[Path] = [
            *BUNDLE_EXTRACTED_MARKERS,
            *HISTORICAL_BUNDLE_EXTRACTED_MARKERS,
        ],
        rerun_files: list[Path] = get_rerun_files(_gas),
        figures_sources: list[Path] = FIGURES_SOURCES,
        outfile: Annotated[Path, Product] = get_methods_figure_file(_gas),
    ) -> None:
        """Generate a gas' methods figure"""
        if gas == "co2":
            generate_co2_figure(outfile, magicc_summary, BUNDLE_DIR)
        elif gas == "ch4":
            generate_ch4_figure(outfile, magicc_summary, BUNDLE_DIR)
        elif gas == "n2o" or gas in N2O_LIKE_GASES:
            generate_n2o_like_figure(gas, outfile, magicc_summary, BUNDLE_DIR)
        elif gas == "cfc12" or gas in CFC12_LIKE_GASES:
            generate_cfc12_like_figure(gas, outfile, BUNDLE_DIR)
        elif gas == "hfc125" or gas in HFC125_LIKE_GASES:
            generate_hfc125_like_figure(gas, outfile, magicc_summary, BUNDLE_DIR)
        else:
            raise NotImplementedError(gas)


for _gas in MAGICC_MISMATCH_GASES:

    @task(id=_gas)
    def task_magicc_mismatch_figure(
        gas: str = _gas,
        magicc_summary: Path = MAGICC_SUMMARY,
        extracted_markers: list[Path] = [
            *BUNDLE_EXTRACTED_MARKERS,
            *HISTORICAL_BUNDLE_EXTRACTED_MARKERS,
        ],
        figures_sources: list[Path] = FIGURES_SOURCES,
        outfile: Annotated[Path, Product] = get_magicc_mismatch_figure_file(_gas),
    ) -> None:
        """Generate a gas' comparison of MAGICC's concentrations with ours"""
        generate_magicc_mismatch_figure(gas, outfile, magicc_summary, BUNDLE_DIR)


for _gas in REGRESSION_GASES:

    @task(id=_gas)
    def task_regressions_figure(
        gas: str = _gas,
        rerun_files: list[Path] = get_rerun_files(_gas),
        figures_sources: list[Path] = FIGURES_SOURCES,
        outfile: Annotated[Path, Product] = get_regressions_figure_file(_gas),
    ) -> None:
        """Generate the figure of the regressions we re-derived for a gas"""
        generate_regressions_figure(gas, outfile, BUNDLE_DIR)


def task_western_et_al_2024_figure(
    extracted_markers: list[Path] = [
        *BUNDLE_EXTRACTED_MARKERS,
        *HISTORICAL_BUNDLE_EXTRACTED_MARKERS,
    ],
    figures_sources: list[Path] = FIGURES_SOURCES,
    outfile: Annotated[Path, Product] = WESTERN_ET_AL_2024_FIGURE,
) -> None:
    """Generate the figure of our extensions of the Western et al. (2024) projections"""
    generate_western_et_al_2024_figure(outfile, BUNDLE_DIR)


for _gas in RESULTS_FIGURE_GASES:

    @task(id=_gas)
    def task_results_figure(
        gas: str = _gas,
        figures_sources: list[Path] = RESULTS_FIGURES_SOURCES,
        outfile: Annotated[Path, Product] = get_results_figure_file(_gas),
    ) -> None:
        """
        Generate a gas' results figure

        The published data it shows comes from ESGF
        (downloaded the first time, then re-used from disk),
        so it isn't tracked as a dependency.
        """
        generate_results_figure(gas, outfile)


for _species in EMISSIONS_FIGURE_SPECIES:

    @task(id=_species)
    def task_emissions_figure(
        species: str = _species,
        extracted_markers: list[Path] = [
            *BUNDLE_EXTRACTED_MARKERS,
            *HISTORICAL_BUNDLE_EXTRACTED_MARKERS,
        ],
        figures_sources: list[Path] = EMISSIONS_FIGURES_SOURCES,
        outfile: Annotated[Path, Product] = get_emissions_figure_file(_species),
    ) -> None:
        """Generate a species' emissions figure"""
        generate_emissions_figure(
            species, outfile, BUNDLE_DIR, cmip6_emissions_file=CMIP6_EMISSIONS_FILE
        )


for _group, _out_file in APPENDIX_FIGURES_FILES.items():

    @task(id=_group)
    def task_appendix_figures(
        group: str = _group,
        texts: list[Path] = [*SECTIONS, *APPENDICES],
        appendix_figures_source: Path = Path(
            local.scenario_ghg_forcing_for_cmip7.appendix_figures.__file__
        ),
        out_file: Annotated[Path, Product] = _out_file,
    ) -> None:
        """Generate the latex for a group of appendix figures"""
        generate_appendix_figures(
            out_file, group=group, figure_tags_in_text=get_tags_in_texts(texts)
        )


def task_author_contributions(
    metadata: Path = METADATA,
    author_contributions_source: Path = Path(
        local.historical_ghg_forcing_for_cmip7.author_contributions.__file__
    ),
    outfile: Annotated[Path, Product] = AUTHOR_CONTRIBUTIONS,
) -> None:
    """Generate the author contributions statement from the metadata file"""
    generate_author_contributions(outfile=outfile, metadata_file=metadata)


if set(APPENDIX_FIGURES_FILES) != set(APPENDIX_FIGURE_GROUPS):
    raise AssertionError
