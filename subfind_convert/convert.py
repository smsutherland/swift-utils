#!/usr/bin/env python

import sys
import os
from yaml import load as yamload, CLoader
import swiftsimio as ssio
from glob import glob
from convert_snap_for_subfind import convert
from make_params import make_param_file

usage = """Usage:
convert.py [-v] [snap num(s)]

Converts one or many swift snapshots into gadget snapshots and sets them up to be used with subfind.
This is meant to be run in the swift directory itself, not the snaps directory.
If you which to convert a single snapshot outside of the context of a swift camels run, use
convert_snap_for_subfind.py"""

def main():
    if len(sys.argv) <= 1:
        print(usage)
        return
    verbose = sys.argv[1] == "-v"

    if not os.path.exists("./snaps/subs"):
        os.mkdir("./snaps/subs")

    try:
        snaps = [int(x) for x in sys.argv[1 + verbose:]]
    except ValueError:
        print(usage)
        return

    for s in snaps:
        convert(s, "./snaps", "./snaps/subs", verbose)

    params = read_from_param_yml()

    params["BoxSize"] = get_boxsize_from_snap("./snaps")

    with open("./snaps/subs/arepo_subfind_params.txt", "w") as f:
        f.write(make_param_file(params))


def read_from_param_yml(fname="params.yml"):
    with open(fname) as f:
        params = yamload(f, CLoader)
    return {
        "Omega_m": params["Cosmology"]["Omega_b"] + params["Cosmology"]["Omega_cdm"],
        "Omega_b": params["Cosmology"]["Omega_b"],
        "H0": params["Cosmology"]["h"],
        "SofteningComovingType1": params["Gravity"]["comoving_DM_softening"],
        "SofteningComoving": params["Gravity"]["comoving_baryon_softening"],
        "SofteningMaxPhysType1": params["Gravity"]["max_physical_DM_softening"],
        "SofteningMaxPhys": params["Gravity"]["max_physical_baryon_softening"],
    }


def get_boxsize_from_snap(snap_dir):
    snaps = glob(snap_dir + "/snapshot_*.hdf5")
    snapname = snaps[0]
    snap = ssio.load(snapname)
    return round( # round to nearest whole number
        (
            snap.metadata.boxsize # in Mpc
            * snap.metadata.cosmology.h # convert to Mpc/h
        ).value.mean() # average over all three axes (they should all be the same anyway)
        * 1000 # convert to kpc/h
    )

if __name__ == "__main__":
    main()
