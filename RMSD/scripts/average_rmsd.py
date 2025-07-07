import os
import argparse
from glob import glob

import pandas as pd
from tqdm import tqdm

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Calculation RMSD averaged over set of trajectories')
    parser.add_argument('--path-to-processing-dirs', required=True)
    parser.add_argument('--output-directory', default=".")
    args = parser.parse_args()

    path_to_processing_dirs = [path_to_processing_dir.lstrip().strip()
                               for path_to_processing_dir in args.path_to_processing_dirs.split(",")]

    time_ns = None
    msd_avg = None

    for path_to_processing_dir in path_to_processing_dirs:
        path_to_rmsd_csv = os.path.join(path_to_processing_dir, "rmsd.csv")

        df_rmsd = pd.read_csv(path_to_rmsd_csv)
        rmsd = df_rmsd[["rmsd_protein", "rmsd_dna", "rmsd_all", "rmsd_dna_inner", "rmsd_dna_outer"]]
        msd = rmsd ** 2

        if time_ns is None:
            time_ns = df_rmsd["time_ns"]

        if msd_avg is None:
            msd_avg = msd
        else:
            msd_avg += msd

    rmsd_avg = (msd_avg / len(path_to_processing_dirs)) ** 0.5
    rmsd_avg["time_ns"] = time_ns

    # Save the rmsd results to a CSV file
    os.makedirs(args.output_directory, exist_ok=True)
    rmsd_avg.to_csv(os.path.join(args.output_directory, "rmsd.csv"), index=False)
