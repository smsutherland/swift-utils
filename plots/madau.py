#!/bin/env python

HELP = """Usage: madau.py [output_fname]
Run in a directory with a SWIFT simulation output and this will create a figure recreating
the SFR density plot from Madau and Dickinson (2014).

Output will be specified by `output_fname`, and will default to "SFR.pdf\""""


def main():
    from astropy.table import QTable
    import astropy.units as u
    import yaml
    import matplotlib.pyplot as plt
    from astropy.visualization import quantity_support
    from sys import argv

    quantity_support()

    sfr = QTable.read(
        "./SFR.txt",
        format="ascii",
        names=[
            "step",
            "Time",
            "a",
            "z",
            "total M_star",
            "SFR (active)",
            "SFR*dt (active)",
            "SFR (total)",
        ],
        units=[
            1.0,
            9.778131e5 * u.Myr,
            1.0,
            1.0,
            1e10 * u.Msun,
            1.02269e-2 * u.Msun / u.yr,
            1e10 * u.Msun,
            1.02269e-2 * u.Msun / u.yr,
        ],
    )
    with open("params.yml", "r") as f:
        params = yaml.safe_load(f)
    h = params["Cosmology"]["h"]
    O_cdm = params["Cosmology"]["Omega_cdm"]
    O_l = params["Cosmology"]["Omega_lambda"]
    O_b = params["Cosmology"]["Omega_b"]
    O_m = O_b + O_cdm
    name = params["MetaData"]["run_name"]

    # TODO: read in the box size from somewhere. It's not in the parameter file
    box_size = 25.0 * u.Mpc / h
    volume = box_size**3
    sfr["SFR density"] = sfr["SFR (total)"] / volume

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.plot(sfr["z"], sfr["SFR density"].to(u.Msun / u.yr / u.Mpc**3), rasterized=True)
    ax.set_yscale("log")
    ax.set_xlim(0, 8)
    ax.set_ylim(10**-2.4, 10**-0.4)
    ax.set_xlabel("redshift (z)")
    ax.set_ylabel(
        rf"SFR density $\left[{u.format.LatexInline.to_string(u.Msun / u.yr / u.Mpc**3)[1:-1]}\right]$"
    )
    ax.set_title(f"{name} SFR density")

    if len(argv) > 1:
        fname = argv[1]
    else:
        fname = "SFR.pdf"

    fig.savefig(fname)


if __name__ == "__main__":
    from sys import argv

    if "-h" in argv or "--help" in argv:
        print(HELP)
        exit(0)
    main()
