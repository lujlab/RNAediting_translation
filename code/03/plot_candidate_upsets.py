#!/usr/bin/env python3
"""Render the original candidate UpSet plots without rerunning site selection."""

import argparse
from pathlib import Path
from rebuild_candidate_corrected_WBG_change2_species_union import draw_plots

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("outdir", type=Path)
    parser.add_argument("n_upset", type=int, nargs="?", default=40)
    args = parser.parse_args()
    draw_plots(args.outdir.resolve(), args.n_upset)
