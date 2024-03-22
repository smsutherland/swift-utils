#!/usr/bin/env python

import sys
import os
from convert_snap_for_subfind import convert

usage = """Usage:
convert.py [-v] [snap num(s)]

Converts one or many swift snapshots into gadget snapshots and sets them up to be used with subfind.
"""

def main():
    if len(sys.argv) <= 1:
        print(usage)
        return
    verbose = sys.argv[1] == "-v"

    if not os.path.exists("./subs"):
        os.mkdir("./subs")

    try:
        snaps = [int(x) for x in sys.argv[1 + verbose:]]
    except ValueError:
        print(usage)
        return

    for s in snaps:
        convert(s, ".", "./subs", verbose)



if __name__ == "__main__":
    main()
