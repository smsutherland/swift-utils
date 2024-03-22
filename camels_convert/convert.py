#!/bin/env python

import sys
from pathlib import Path
import cosmoconvert.gadget2swift as convert
from param_file import make_param_file
import shutil

usage = """Usage:
    ./convert.py [CAMELS sim directory] [target SWIFT IC] [name]"""


def main():
    Omega_b = 0.049
    argv = sys.argv[1:]

    if len(argv) != 3:
        print(usage)
        sys.exit(1)

    [sim_dir, target_dir, run_name] = argv
    sim_dir = Path(sim_dir)
    target_dir = Path(target_dir)

    convert_ic(sim_dir / "ICs", target_dir / "ic.hdf5", omega_baryon=Omega_b)

    params = read_params(sim_dir)
    # All SIMBA runs have Omega_b=0.049 and h=0.6711
    write_param_file(target_dir, Omega_b=Omega_b, run_name=run_name, **params)

    copy_misc_files(target_dir)
    write_cosmo_params(target_dir, params)

def convert_ic(ic_dir: Path, target: Path, omega_baryon):
    convert.gadget2swift(str(ic_dir / "ics"), omega_baryon, out=target)

def write_param_file(target_dir: Path, **kwargs):
    with open(target_dir/"params.yml", "w") as f:
        f.write(make_param_file(**kwargs))

def read_params(sim_dir: Path) -> dict:
    with open(sim_dir / "CosmoAstro_params.txt") as f:
        params = [float(x) for x in f.read().split()]
    names = [
        "Omega_m",
        "sigma_8",
        "A_SN1",
        "A_AGN1",
        "A_SN2",
        "A_AGN2",
    ]
    return dict(zip(names, params))

def copy_misc_files(target_dir: Path):
    dir = Path(__file__).parent / "copy_files"
    for fname in dir.iterdir():
        if fname.is_file():
            shutil.copy(fname, target_dir)
        else:
            shutil.copytree(fname, target_dir / fname.name)

    dir = Path(__file__).parent / "link_files"
    for fname in dir.iterdir():
        (target_dir / fname.name).symlink_to(fname)

def write_cosmo_params(target_dir: Path, params: dict):
    with open(target_dir / "CosmoAstro_params.txt", "w") as f:
        f.write(
            f"{params['Omega_m']:.5f} "
            f"{params['sigma_8']:.5f} "
            f"{params['A_SN1']:.5f} "
            f"{params['A_AGN1']:.5f} "
            f"{params['A_SN2']:.5f} "
            f"{params['A_AGN2']:.5f}\n"
        )


if __name__ == "__main__":
    main()
