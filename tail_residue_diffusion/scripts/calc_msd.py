import argparse
import numpy as np
import os
import pandas as pd

from glob import glob

from process_utils.calc import calc_mean_square_displacement

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate mean square displacement of center of mass (msd)')

    parser.add_argument('--atom-coords-dir')
    parser.add_argument('--output-directory')
    parser.add_argument('--lag-spacing', choices=["linear", "log"], default="log")
    parser.add_argument('--n-lag-points', default=1000, type=str)
    parser.add_argument('--trj-length', default=10, type=int)
    args = parser.parse_args()

    atom_coords_csvs = glob(os.path.join(args.atom_coords_dir, "*csv"))
    atom_coords_csvs.sort()

    for atom_coords_csv in atom_coords_csvs:

        atom_coords = pd.read_csv(atom_coords_csv)
        time_ns = atom_coords["time_ns"]
        xyz = atom_coords[["x", "y", "z"]].values

        # if args.lag_spacing == "linear":
        #     lag_index = np.unique(np.linspace(1, len(time_ns), args.n_lag_points, endpoint=False).astype(int))
        # else:
        #     lag_index = np.unique(np.logspace(0, np.log10(len(time_ns) - 1), args.n_lag_points).astype(int))

        lag, msd = calc_mean_square_displacement(time_ns, xyz, lag_index=range(1, len(time_ns)))

        os.makedirs(args.output_directory, exist_ok=True)
        pd.DataFrame({
            "time_ns": [0, *lag],
            "msd": [0, *msd],
        }).to_csv(os.path.join(args.output_directory, os.path.basename(atom_coords_csv)), index=False)
