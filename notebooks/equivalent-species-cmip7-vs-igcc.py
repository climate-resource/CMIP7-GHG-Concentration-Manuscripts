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
#   title: Equivalent species, CMIP7 vs. IGCC
# ---

# %% [markdown]
# # Equivalent species, CMIP7 vs. IGCC
#
# Here we show why our (CMIP7) equivalent species differ from those of
# the Indicators of Global Climate Change (IGCC) assessment
# (Forster et al., 2026, https://doi.org/10.5194/essd-18-3889-2026)
# and check the statements about this in `results.tex`.
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
# IGCC publishes CFC-12 and HFC-134a equivalents, but no CFC-11 equivalent,
# so there is no CFC-11 equivalent comparison here.
#
# In short, the answer is that
#
# - for CFC-12 equivalent, the radiative efficiencies drive most of the difference.
#   IGCC uses those of Hodnebrog et al. (2020), we use those of AR6.
#   These are almost the same, except that AR6's are higher for CFC-11 and CFC-12,
#   so in AR6 every other gas has less weight relative to CFC-12.
#   IGCC also includes eight gases which we don't, which drives most of the rest.
# - for HFC-134a equivalent, the gases included drive most of the difference.
#   IGCC includes only the HFCs, we also include the PFCs, SF6, NF3 and SO2F2.

# %% [markdown]
# ## Imports

# %%
import pandas as pd

from local.historical_ghg_forcing_for_cmip7.equivalent_species import (
    LatexClaim,
    assert_claims_hold,
    check_reproduction,
    decompose_difference,
    get_contributions_by_gas,
    get_radiative_effect,
    load_cmip7,
    load_igcc,
    plot_decomposition,
    summarise_claims,
)

# %%
pd.set_option("display.max_columns", 20)
pd.set_option("display.width", 200)
pd.set_option("display.max_colwidth", 120)

# %%
EQUIVALENT_SPECIES = ("cfc12eq", "hfc134aeq")
"""The equivalent species both datasets have"""

# %% [markdown]
# ## Load data

# %%
cmip7 = load_cmip7()
igcc = load_igcc()

# %% [markdown]
# ## Check that we have the definitions right
#
# If the definitions reproduce each dataset's published equivalent species,
# we know we have copied them correctly.

# %%
reproduction = pd.DataFrame(
    {
        (dataset.label, equivalent_species): {
            "max abs residual": check_reproduction(equivalent_species, dataset)
            .abs()
            .max()
        }
        for dataset in (cmip7, igcc)
        for equivalent_species in EQUIVALENT_SPECIES
    }
).T
reproduction

# %%
REPRODUCTION_TOLERANCE = 1e-9
for dataset in (cmip7, igcc):
    for equivalent_species in EQUIVALENT_SPECIES:
        residual = check_reproduction(equivalent_species, dataset).abs().max()
        assert residual < REPRODUCTION_TOLERANCE, (
            dataset.label,
            equivalent_species,
            residual,
        )

# %% [markdown]
# ## Gases included

# %%
for equivalent_species in EQUIVALENT_SPECIES:
    ours = set(cmip7.definition.components[equivalent_species])
    theirs = set(igcc.definition.components[equivalent_species])
    print(f"{equivalent_species}")
    print(f"- only in CMIP7: {sorted(ours - theirs)}")
    print(f"- only in IGCC: {sorted(theirs - ours)}")

# %% [markdown]
# ## Radiative efficiencies
#
# What matters for an equivalent species is each gas' radiative efficiency
# relative to the reference gas' (its "weight").
# AR6's radiative efficiencies of CFC-11 and CFC-12 are around 12% higher
# than those of Hodnebrog et al. (2020), while the others are almost the same.

# %%
radiative_efficiencies = pd.DataFrame(
    {
        "CMIP7 (AR6) [W / m^2 / ppb]": cmip7.definition.radiative_efficiencies,
        "IGCC (Hodnebrog et al., 2020) [W / m^2 / ppb]": igcc.definition.radiative_efficiencies,
    }
).dropna()
radiative_efficiencies["ratio"] = (
    radiative_efficiencies["CMIP7 (AR6) [W / m^2 / ppb]"]
    / radiative_efficiencies["IGCC (Hodnebrog et al., 2020) [W / m^2 / ppb]"]
)
radiative_efficiencies.sort_values("ratio")

# %% [markdown]
# ## Differences
#
# We look at each equivalent species in turn.
# All differences are CMIP7 minus IGCC.
# IGCC starts in 1850 (its 1750 value is a pre-industrial reference)
# and our dataset ends before IGCC's,
# so the comparison covers the years in between.

# %%
decompositions = {
    equivalent_species: decompose_difference(equivalent_species, cmip7, igcc)
    for equivalent_species in EQUIVALENT_SPECIES
}
last_year = min(d.index.max() for d in decompositions.values())
last_year

# %%
SELECTED_YEARS = [1850, 1950, 1980, 1990, 2000, 2010, 2014, 2020, last_year]


def show_equivalent_species(equivalent_species: str) -> None:
    """Show everything about the difference in one equivalent species"""
    decomposition = decompositions[equivalent_species]
    plot_decomposition(equivalent_species, decomposition, cmip7, igcc)

    display(decomposition.loc[SELECTED_YEARS].round(2))  # noqa: F821

    print(f"Contributions to the difference by gas in {last_year}:")
    display(  # noqa: F821
        get_contributions_by_gas(equivalent_species, cmip7, igcc, last_year)
        .head(12)
        .round(2)
    )


# %% [markdown]
# ### CFC-12 equivalent
#
# Most of the difference is from the radiative efficiencies,
# mainly through the weight given to HCFC-22, CFC-113 and CCl4.
# Most of the rest is from the gases only IGCC includes,
# the largest being CFC-13.

# %%
show_equivalent_species("cfc12eq")

# %% [markdown]
# ### HFC-134a equivalent
#
# Almost all of the difference is from the gases only we include,
# mainly SF6 and CF4.
# Most of the rest is from differences in the concentrations of HFC-134a and HFC-32.

# %%
show_equivalent_species("hfc134aeq")

# %% [markdown]
# ## Check the statements in `results.tex`
#
# Each statement about the IGCC comparison is written out below
# with the range of values for which it holds.
# Statements about how big the difference is "now"
# are checked against the last year both datasets cover.
# If any doesn't hold, the last cell fails:
# update `results.tex` (and the statement here) to match the data.

# %%
LATEX_FILE = "manuscripts/historical-ghg-forcing-for-cmip7/results.tex"

cfc12eq = decompositions["cfc12eq"]
cfc12eq_now = cfc12eq.loc[last_year, "total difference"]
cfc12eq_max = cfc12eq["total difference"].abs().max()

hfc134aeq = decompositions["hfc134aeq"]
hfc134aeq_now = hfc134aeq.loc[last_year, "total difference"]
hfc134aeq_max = hfc134aeq["total difference"].abs().max()

claims = [
    # TODO: consider whether this can be done more easily
    LatexClaim(
        section="sssec:results-equivalence-datasets-cfc12eq",
        stated="Our dataset is around 40~ppt lower than the IGCC dataset",
        statistic=f"CMIP7 - IGCC in {last_year} [ppt]",
        calculated=cfc12eq_now,
        lower=-45.0,
        upper=-35.0,
    ),
    LatexClaim(
        section="sssec:results-equivalence-datasets-cfc12eq",
        stated="(~0.015~W/m2)",
        statistic="approx. radiative effect of the above [W / m^2]",
        calculated=get_radiative_effect("cfc12eq", abs(cfc12eq_now)),
        lower=0.010,
        upper=0.020,
    ),
    LatexClaim(
        section="sssec:results-equivalence-datasets-cfc12eq",
        stated=(
            "the difference in radiative efficiencies is the dominant one, "
            "contributing approximately 90\\% of the difference"
        ),
        statistic=f"radiative efficiency share of the difference in {last_year}",
        calculated=cfc12eq.loc[last_year, "difference in radiative efficiency"]
        / cfc12eq_now,
        lower=0.85,
        upper=0.95,
    ),
    LatexClaim(
        section="sssec:results-equivalence-datasets-cfc12eq",
        stated="in approximate radiative forcing terms the difference is about 0.015~W/m2",
        statistic="approx. radiative effect of max abs CMIP7 - IGCC [W / m^2]",
        calculated=get_radiative_effect("cfc12eq", cfc12eq_max),
        lower=0.01,
        upper=0.02,
    ),
    LatexClaim(
        section="sssec:results-equivalence-datasets-hfc134aeq",
        stated="Our dataset is around 125~ppt higher than the IGCC dataset",
        statistic=f"CMIP7 - IGCC in {last_year} [ppt]",
        calculated=hfc134aeq_now,
        lower=115.0,
        upper=135.0,
    ),
    LatexClaim(
        section="sssec:results-equivalence-datasets-hfc134aeq",
        stated="(~0.02~W/m2)",
        statistic="approx. radiative effect of the above [W / m^2]",
        calculated=get_radiative_effect("hfc134aeq", hfc134aeq_now),
        lower=0.015,
        upper=0.025,
    ),
    LatexClaim(
        section="sssec:results-equivalence-datasets-hfc134aeq",
        stated=(
            "the IGCC's exclusion of any gas which is not an HFC [...] is the dominant one, "
            "contributing approximately 90\\% of the difference"
        ),
        statistic=f"included gases share of the difference in {last_year}",
        calculated=hfc134aeq.loc[last_year, "difference in included gases"]
        / hfc134aeq_now,
        lower=0.85,
        upper=0.95,
    ),
    LatexClaim(
        section="sssec:results-equivalence-datasets-hfc134aeq",
        stated="in approximate radiative forcing terms the difference is less than 0.05~W/m2",
        statistic="approx. radiative effect of max abs CMIP7 - IGCC [W / m^2]",
        calculated=get_radiative_effect("hfc134aeq", hfc134aeq_max),
        upper=0.05,
    ),
]
summarise_claims(claims)

# %%
assert_claims_hold(claims, LATEX_FILE)
