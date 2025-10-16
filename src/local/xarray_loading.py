"""
Loading of xarray objects
"""

from __future__ import annotations

from collections.abc import Iterable
from functools import partial
from pathlib import Path
from typing import Callable

import ncdata
import xarray as xr
import xarray.backends.api
from ncdata.netcdf4 import from_nc4
from ncdata.xarray import to_xarray

from local.esgf.models import ESGFDataset


def load_xarray_from_esgf_dataset(  # noqa: PLR0912, PLR0913
    esgf_dataset: ESGFDataset,
    pre_to_xarray: Callable[
        [list[ncdata.NcData], ESGFDataset, list[Path]], list[ncdata.NcData]
    ]
    | None = None,
    concat_dim: str = "time",
    post_to_xarray: Callable[[xr.Dataset, ESGFDataset, list[Path]], xr.Dataset]
    | None = None,
    map_variable_name: bool = True,
    add_auxilliary_sources: bool | tuple[str, ...] = False,
    map_auxilliary_sources_names: bool = True,
    add_attributes_from_metadata: Iterable[str] | None = None,
    load_auxilliary_array_from_esgf_dataset: Callable[[ESGFDataset], xr.Dataset]
    | None = None,
    merge_auxilliary_ds: Callable[[list[xr.Dataset]], xr.Dataset] | None = None,
    create_default_indexes: bool = True,
) -> xr.Dataset:
    """
    Load an [xr.Dataset][xarray.Dataset] from a [ESGFDataset][local.models.esgf.]

    Parameters
    ----------
    esgf_dataset
        [ESGFDataset][local.models.esgf.] from which to load

    pre_to_xarray
        Processing to apply before loading with [xarray][]

        If supplied, this must be a callable which receives [ncdata][] objects,
        the `esgf_dataset` and the associated filepaths
        and returns [ncdata][] objects.
        It can do whatever manipulations are needed.

    concat_dim
        Dimension along which to concatenate the loaded files

        Only used if there is more than one file associated with `esgf_dataset`

    post_to_xarray
        Processing to apply after converting to [xr.Dataset][xarray.Dataset]

        If supplied,
        this must be a callable which receives the [xr.Dataset][xarray.Dataset] object,
        the `esgf_dataset` and the associated filepaths
        and returns an [xr.Dataset][xarray.Dataset] object.
        It can do whatever manipulations are needed.

    map_variable_name
        If `True`, we map from the original variable name in the file
        to the variable name as given in `esgf_dataset`

    add_auxilliary_sources
        If `True`, add all auxilliary sources to the loaded data

        If a tuple of strings, only add given auxilliary sources to the loaded data.

        If `False`, don't add auxilliary sources to the loaded data

    map_auxilliary_sources_names
        If `True`, we map from the original variable name in the auxilliary files
        to the corrected variable names

    add_attributes_from_metadata
        If supplied, attributes to set on the result from values in `esgf_dataset`

    load_auxilliary_array_from_esgf_dataset
        Callable to use to load auxilliary data

        If not supplied, we use [load_xarray_from_esgf_dataset][(m).]

    merge_auxilliary_ds
        Callable to use to merge auxilliary data with the 'main' data

        If not supplied, we use `partial(xr.merge, compat="broadcast_equals")`

    create_default_indexes
        Create default indexes on the result

    Returns
    -------
    :
        Loaded [xr.Dataset][xarray.Dataset]
    """
    if load_auxilliary_array_from_esgf_dataset is None:
        load_auxilliary_array_from_esgf_dataset = partial(
            load_xarray_from_esgf_dataset,
            map_variable_name=map_auxilliary_sources_names,
        )

    if merge_auxilliary_ds is None:
        merge_auxilliary_ds = partial(xr.merge, compat="broadcast_equals")

    fps = [Path(f.esgf_file_local.path) for f in esgf_dataset.esgf_files]
    ncdatas = [from_nc4(fp) for fp in fps]

    if pre_to_xarray is not None:
        ncdatas = pre_to_xarray(ncdatas, esgf_dataset, fps)

    time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)
    to_xr = partial(to_xarray, decode_times=time_coder)
    if len(ncdatas) > 1:
        ds = xr.concat(
            [to_xarray(nc, decode_times=time_coder) for nc in ncdatas],
            dim=concat_dim,
            data_vars=None,
        )

    else:
        ds = to_xr(ncdatas[0])

    if create_default_indexes:
        # Ensure that the dataset has the right indexes too
        # so that later operations (e.g. selecting) work.
        ds = xarray.backends.api._maybe_create_default_indexes(ds)
        # Note:
        # The above is how xarray works.
        # We have to add this in a custom way here
        # because ncdata's `to_xarray` function uses a route
        # that doesn't support this argument
        # (maybe missing implementation here
        # because the `set_indexes` argument isn't used:
        # https://github.com/pydata/xarray/blob/3572f4e70f2b12ef9935c1f8c3c1b74045d2a092/xarray/backends/store.py#L36C9-L36C20
        # ?)

    if post_to_xarray is not None:
        ds = post_to_xarray(ds, esgf_dataset, fps)

    variable_mappings = {}
    if map_variable_name:
        variable_mappings[esgf_dataset.esgf_raw_metadata.variable] = (
            esgf_dataset.variable
        )

    if add_auxilliary_sources:
        raise NotImplementedError
        for aux_source in esgf_dataset.auxilliary_sources:
            if aux_source.internal:
                aux_variable = aux_source.internal.variable_name_data

            else:
                aux_variable = aux_source.external.variable_root_dd

            if (
                isinstance(add_auxilliary_sources, bool) and add_auxilliary_sources
            ) or (
                isinstance(add_auxilliary_sources, tuple)
                and aux_variable in add_auxilliary_sources
            ):
                if aux_source.internal:
                    if map_auxilliary_sources_names:
                        variable_mappings[aux_source.internal.variable_name_data] = (
                            aux_source.internal.variable_root_dd
                        )

                else:
                    auxilliary_ds = load_auxilliary_array_from_esgf_dataset(
                        esgf_dataset=aux_source.external,
                    )

                    ds = merge_auxilliary_ds([ds, auxilliary_ds])

    ds = ds.rename_vars(variable_mappings)

    if add_attributes_from_metadata is not None:
        for attr_to_add in add_attributes_from_metadata:
            ds = ds.assign_attrs({attr_to_add: getattr(esgf_dataset, attr_to_add)})

    return ds
