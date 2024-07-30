#!/bin/env python

import sys

from matplotlib.colors import LogNorm
import h5py
import matplotlib.pyplot as plt
import numpy as np
import swiftsimio as ssio
import unyt as u
from swiftsimio.visualisation.projection import project_gas, project_pixel_grid
from swiftsimio.visualisation.smoothing_length_generation import (
    generate_smoothing_lengths,
)


def main():
    [snap, fof, nth] = sys.argv[1:]
    nth = int(nth)
    data = ssio.load(snap)
    with h5py.File(fof) as fof:
        masses = fof["Group/GroupMass"][:]
        group_idx = np.argsort(masses)[-nth]
        group_idx = nth
        center = fof["Group/GroupPos"][group_idx] / data.metadata.cosmology.h
        radius = fof["Group/Group_R_Crit200"][group_idx] / data.metadata.cosmology.h

    r_factor = 1.5
    e = radius * r_factor
    minimum = center - e
    maximum = center + e
    region = np.array(
        [
            minimum[0],
            maximum[0],
            minimum[1],
            maximum[1],
            minimum[2],
            maximum[2],
        ]
    )
    region = ssio.objects.cosmo_array(region, u.kpc, cosmo_factor=data.metadata.a)

    gas_mass = project_gas(
        data,
        resolution=1024,
        project="masses",
        parallel=True,
        periodic=True,
        region=region,
    )

    data.dark_matter.smoothing_length = generate_smoothing_lengths(
        data.dark_matter.coordinates,
        data.metadata.boxsize,
        kernel_gamma=1.8,
        neighbours=57,
        speedup_fac=2,
        dimension=3,
    )
    dm_mass = (
        project_pixel_grid(
            data=data.dark_matter,
            boxsize=data.metadata.boxsize,
            resolution=1024,
            project="masses",
            parallel=True,
            region=region,
            periodic=True,
        )
        * gas_mass.units
    )

    clip_factor = 10
    gas_mass_clipped = np.clip(gas_mass.value, 0, clip_factor * np.mean(gas_mass.value))
    dm_mass_clipped = np.clip(dm_mass.value, 0, clip_factor * np.mean(dm_mass.value))
    rgb = np.stack(
        [
            # LogNorm()(mass_grid),
            # LogNorm()(dm_mass),
            gas_mass_clipped / np.max(gas_mass_clipped),
            dm_mass_clipped / np.max(dm_mass_clipped),
            np.zeros(dm_mass.shape),
        ],
        axis=-1,
    )
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.imshow(rgb, extent=[-e, e, -e, e])
    ax.add_patch(plt.Circle((0, 0), radius, color="w", fill=False))
    ax.set_xlabel("x [kpc/h]")
    ax.set_ylabel("y [kpc/h]")
    fig.tight_layout()
    fig.savefig(f"halo_{nth}.pdf")


if __name__ == "__main__":
    main()
