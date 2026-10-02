# ---
# jupyter:
#   authors:
#   - name: Zebedee Nicholls
#   jupytext:
#     notebook_metadata_filter: title,authors
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.17.3
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
#   title: Equivalent species, CMIP7 vs. CMIP6
# ---

# %% [markdown]
# # Equivalent species, CMIP7 vs. CMIP6
#
# Here we show why our (CMIP7) equivalent species differ from CMIP6's
# and check the statements about this in the historical manuscript.
#
# An equivalent species is a radiative efficiency weighted sum
# of the concentrations of a group of gases.
# Two datasets' equivalent species can therefore differ because of
#
# - a difference in the gases included
# - a difference in the radiative efficiencies used
# - other differences, i.e. differences in the concentrations of the individual gases
#
# For how we split the difference into these parts,
# see `local.historical_ghg_forcing_for_cmip7.equivalent_species`.
#
# In short, the answer is that CMIP6 includes the same gases as us,
# but used radiative efficiencies from the IPCC's Fifth Assessment Report (AR5)
# whereas we use those from the Sixth (AR6).
# This, not the concentrations, is what drives the differences
# in CFC-11 and CFC-12 equivalent.

# %% [markdown]
# ## Imports

# %%
import pandas as pd

from local.historical_ghg_forcing_for_cmip7.equivalent_species import (
    EQUIVALENT_SPECIES,
    UNITS,
    LatexClaim,
    assert_claims_hold,
    check_reproduction,
    decompose_difference,
    get_contributions_by_gas,
    get_radiative_effect,
    load_cmip6,
    load_cmip7,
    plot_decomposition,
    summarise_claims,
)

# %%
pd.set_option("display.max_columns", 20)
pd.set_option("display.width", 200)
pd.set_option("display.max_colwidth", 120)

# %% [markdown]
# ## Load data

# %%
cmip7 = load_cmip7()
cmip6 = load_cmip6()

# %% [markdown]
# ## Check that we have the definitions right
#
# If the definitions reproduce each dataset's published equivalent species,
# we know we have copied them correctly.
# The CMIP6 data is stored as 32-bit floats, hence the looser tolerance.

# %%
reproduction = pd.DataFrame(
    {
        (dataset.label, equivalent_species): {
            "max abs residual": check_reproduction(equivalent_species, dataset)
            .abs()
            .max()
        }
        for dataset in (cmip7, cmip6)
        for equivalent_species in EQUIVALENT_SPECIES
    }
).T
reproduction

# %%
for dataset, tolerance in ((cmip7, 1e-9), (cmip6, 1e-3)):
    for equivalent_species in EQUIVALENT_SPECIES:
        residual = check_reproduction(equivalent_species, dataset).abs().max()
        assert residual < tolerance, (dataset.label, equivalent_species, residual)

# %% [markdown]
# This also confirms that the CMIP6 equivalent species include each gas once
# and only once: double counting or missing a gas would show up here.
#
# CMIP6 includes the same gases as us in each equivalent species.

# %%
for equivalent_species in EQUIVALENT_SPECIES:
    assert set(cmip7.definition.components[equivalent_species]) == set(
        cmip6.definition.components[equivalent_species]
    ), equivalent_species

# %% [markdown]
# ## Radiative efficiencies
#
# What matters for an equivalent species is each gas' radiative efficiency
# relative to the reference gas' (its "weight").
# AR6's radiative efficiencies of CFC-11 and CFC-12 are higher than AR5's
# by more than most other gases',
# so in AR6 every other gas has less weight relative to CFC-11 and CFC-12.

# %%
radiative_efficiencies = pd.DataFrame(
    {
        "CMIP7 (AR6) [W / m^2 / ppb]": cmip7.definition.radiative_efficiencies,
        "CMIP6 (AR5) [W / m^2 / ppb]": cmip6.definition.radiative_efficiencies,
    }
).dropna()
radiative_efficiencies["ratio"] = (
    radiative_efficiencies["CMIP7 (AR6) [W / m^2 / ppb]"]
    / radiative_efficiencies["CMIP6 (AR5) [W / m^2 / ppb]"]
)
radiative_efficiencies.sort_values("ratio")

# %% [markdown]
# ## Differences
#
# We look at each equivalent species in turn.
# All differences are CMIP7 minus CMIP6.

# %%
decompositions = {
    equivalent_species: decompose_difference(equivalent_species, cmip7, cmip6)
    for equivalent_species in EQUIVALENT_SPECIES
}

# %%
SELECTED_YEARS = [1, 1750, 1850, 1950, 1980, 1990, 2000, 2010, 2014]


def show_equivalent_species(equivalent_species: str) -> None:
    """Show everything about the difference in one equivalent species"""
    decomposition = decompositions[equivalent_species]
    plot_decomposition(equivalent_species, decomposition, cmip7, cmip6)

    display(decomposition.loc[SELECTED_YEARS].round(2))  # noqa: F821

    max_deviation_year = decomposition["total difference"].abs().idxmax()
    print(
        f"Maximum deviation: {decomposition.loc[max_deviation_year, 'total difference']:.2f} {UNITS} "
        f"in {max_deviation_year}"
    )
    print(
        "Contributions to the difference by gas "
        f"in the year of maximum deviation ({max_deviation_year}):"
    )
    display(  # noqa: F821
        get_contributions_by_gas(equivalent_species, cmip7, cmip6, max_deviation_year)
        .head(10)
        .round(2)
    )


# %% [markdown]
# ### CFC-12 equivalent
#
# The difference is mostly due to the radiative efficiencies
# (at least three quarters of it in every year).
# Before the industrial era, it is entirely due to them:
# the concentrations are practically the same,
# but AR6's much lower radiative efficiency of CH3Cl relative to CFC-12
# almost halves the CFC-12 equivalent.

# %%
show_equivalent_species("cfc12eq")

# %% [markdown]
# ### CFC-11 equivalent
#
# The same story as CFC-12 equivalent.

# %%
show_equivalent_species("cfc11eq")

# %% [markdown]
# ### HFC-134a equivalent
#
# The differences are much smaller.
# Before the industrial era, the difference is entirely due to the radiative efficiencies.
# Since then, differences in concentrations matter about as much.

# %%
show_equivalent_species("hfc134aeq")

# %% [markdown]
# ## Check the statements in `results.tex`
#
# Each statement about the CMIP6 comparison is written out below
# with the range of values for which it holds.
# If any doesn't hold, the last cell fails:
# update `results.tex` (and the statement here) to match the data.

# %%
LATEX_FILE = "manuscripts/historical-ghg-forcing-for-cmip7/results.tex"

cfc12eq = decompositions["cfc12eq"]
cfc12eq_since_1995 = cfc12eq.loc[1995:, "total difference"].abs().max()
cfc12eq_max_year = cfc12eq["total difference"].abs().idxmax()
# The difference from CMIP6 recalculated with AR6 radiative efficiencies
# is everything other than the radiative efficiency part
cfc12eq_recalculated = (
    (cfc12eq["total difference"] - cfc12eq["difference in radiative efficiency"])
    .abs()
    .max()
)

hfc134aeq = decompositions["hfc134aeq"]
hfc134aeq_max_year = hfc134aeq["total difference"].abs().idxmax()
hfc134aeq_max = hfc134aeq["total difference"].abs().max()

claims = [
    LatexClaim(
        section="sssec:results-equivalence-datasets-cfc12eq",
        stated="the maximum deviation is 40~ppt from 1995 onwards",
        statistic="max abs CMIP7 - CMIP6, 1995 onwards [ppt]",
        calculated=cfc12eq_since_1995,
        lower=35.0,
        upper=45.0,
    ),
    LatexClaim(
        section="sssec:results-equivalence-datasets-cfc12eq",
        stated="(~0.015~W/m2)",
        statistic="approx. radiative effect of the above [W / m^2]",
        calculated=get_radiative_effect("cfc12eq", cfc12eq_since_1995),
        lower=0.0125,
        upper=0.0175,
    ),
    LatexClaim(
        section="sssec:results-equivalence-datasets-cfc12eq",
        stated="The changes from CMIP6 are driven by the use of updated radiative efficiencies",
        statistic="radiative efficiency share of the difference in the year of maximum deviation",
        calculated=cfc12eq.loc[cfc12eq_max_year, "difference in radiative efficiency"]
        / cfc12eq.loc[cfc12eq_max_year, "total difference"],
        lower=0.5,
    ),
    LatexClaim(
        section="sssec:results-equivalence-datasets-cfc12eq",
        stated=(
            "If the CMIP6 CFC12-eq is recalculated using AR6 radiative efficiencies, "
            "the difference from our dataset reduces to below 8~ppt"
        ),
        statistic="max abs CMIP7 - CMIP6 recalculated with AR6 radiative efficiencies [ppt]",
        calculated=cfc12eq_recalculated,
        upper=8.0,
    ),
    LatexClaim(
        section="sssec:results-equivalence-datasets-cfc12eq",
        stated="(~0.003~W/m2)",
        statistic="approx. radiative effect of the above [W / m^2]",
        calculated=get_radiative_effect("cfc12eq", cfc12eq_recalculated),
        lower=0.0025,
        upper=0.0035,
    ),
    LatexClaim(
        section="sssec:results-equivalence-datasets-hfc134aeq",
        stated="the maximum deviation is 2.5~ppt",
        statistic="max abs CMIP7 - CMIP6 [ppt]",
        calculated=hfc134aeq_max,
        lower=2.0,
        upper=3.0,
    ),
    LatexClaim(
        section="sssec:results-equivalence-datasets-hfc134aeq",
        stated="in 1980",
        statistic="year of maximum deviation",
        calculated=hfc134aeq_max_year,
        lower=1975,
        upper=1985,
    ),
    LatexClaim(
        section="sssec:results-equivalence-datasets-hfc134aeq",
        stated="(~0.0005~W/m2)",
        statistic="approx. radiative effect of the maximum deviation [W / m^2]",
        calculated=get_radiative_effect("hfc134aeq", hfc134aeq_max),
        lower=0.0004,
        upper=0.0006,
    ),
    LatexClaim(
        section="sssec:results-equivalence-datasets-hfc134aeq",
        stated="This is driven by updates to the input data used to create the dataset",
        statistic='"other" (i.e. concentration) share of the difference in the year of maximum deviation',
        calculated=hfc134aeq.loc[hfc134aeq_max_year, "other"]
        / hfc134aeq.loc[hfc134aeq_max_year, "total difference"],
        lower=0.5,
    ),
]
summarise_claims(claims)

# %%
assert_claims_hold(claims, LATEX_FILE)
