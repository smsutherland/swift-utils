#!/bin/env python

import yt
import sys

help = """USAGE:
    slice.py [snapshot name] [output name]
"""


def main():
    if len(sys.argv) != 3:
        print(help)
        return

    [_, snap_name, out_name] = sys.argv
    make_slice(snap_name, out_name)


def make_slice(snap_name: str, out_name: str):
    ds = yt.load(snap_name, hint="swift")
    tprj = yt.SlicePlot(
        ds,
        "x",
        ("gas", "temperature"),
        # width=(150, "kpc"),
        fontsize=15,
    )
    # tprj.set_zlim(("gas", "temperature"), 1e2, 1e9)
    tprj.set_cmap(("gas", "temperature"), "dusk")
    tprj.set_colorbar_label("temperature", "K")
    tprj.set_axes_unit("kpc")
    tprj.save(out_name, "pdf")


if __name__ == "__main__":
    main()
