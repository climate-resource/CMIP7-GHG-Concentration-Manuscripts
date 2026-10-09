"""
Functionality related to generating the historical GHG manuscript
"""

from .appendix_figures import generate_appendix_figures, get_figure_tags_in_text
from .c4f10_like_methods_figure import (
    C4F10_LIKE_GASES,
    generate_c4f10_like_methods_figure,
)
from .c8f18_methods_figure import generate_c8f18_methods_figure
from .cfc12_like_methods_figure import (
    CFC12_LIKE_GASES,
    generate_cfc12_like_methods_figure,
)
from .cfc12_like_tables import (
    generate_cfc12_like_obs_network_sources_list,
    generate_cfc12_like_per_gas_table,
)
from .ch4_methods_figure import generate_ch4_methods_figure
from .co2_methods_figure import generate_co2_methods_figure
from .input_data_references_tables import generate_input_data_references_tables
from .n2o_methods_figure import generate_n2o_methods_figure
from .results_figure import generate_results_figure_for_gas

__all__ = [
    "C4F10_LIKE_GASES",
    "CFC12_LIKE_GASES",
    "generate_appendix_figures",
    "generate_c4f10_like_methods_figure",
    "generate_c8f18_methods_figure",
    "generate_cfc12_like_methods_figure",
    "generate_cfc12_like_obs_network_sources_list",
    "generate_cfc12_like_per_gas_table",
    "generate_ch4_methods_figure",
    "generate_co2_methods_figure",
    "generate_input_data_references_tables",
    "generate_n2o_methods_figure",
    "generate_results_figure_for_gas",
    "get_figure_tags_in_text",
]
