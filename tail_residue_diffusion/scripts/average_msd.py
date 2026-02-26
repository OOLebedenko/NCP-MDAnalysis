import argparse
import numpy as np
import os
import pandas as pd

from glob import glob
from tqdm import tqdm


def average_across_trajectories(arrays, dt_ns=0.01) -> tuple:
    lengths = np.array([len(arr) for arr in arrays])
    max_lag = np.max(lengths)

    # Create weight matrix with padding for shorter trajectories
    k = np.arange(max_lag)
    weights = np.maximum(0, lengths[:, np.newaxis] - k)  # shape: (n_trajectories, max_lag)

    # Create values matrix with zeros where trajectory is shorter
    values = np.zeros((len(arrays), max_lag))
    for i, (arr, N) in enumerate(zip(arrays, lengths)):
        values[i, :N] = arr[:N]

    # Weighted sum (weights are zero for missing data)
    weighted_sum = np.sum(values * weights, axis=0)
    total_weights = np.sum(weights, axis=0)

    # add time data
    time_ns = np.arange(max_lag) * dt_ns

    # Avoid division by zero
    return time_ns, np.divide(weighted_sum, total_weights, where=total_weights > 0)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Average msd over trajectories of the different length')
    parser.add_argument('--path-to-processing-dirs', required=True)
    parser.add_argument('--output-directory', default=".")
    args = parser.parse_args()

    path_to_processing_dirs = [path_to_processing_dir.lstrip().strip()
                               for path_to_processing_dir in args.path_to_processing_dirs.split(",")]

    file_basenames = [os.path.basename(acorr_csv).split("-", 1)[-1] for acorr_csv in
                      glob(os.path.join(path_to_processing_dirs[0], "data", "msd", "B*.csv"))]
    file_basenames.sort()

    tail_copies = ["B", "F"]
    dt_ns = None
    for file_basename in tqdm(file_basenames, desc="average msd"):
        acorr_arrays = []
        for path_to_processing_dir in path_to_processing_dirs:
            for tail_copy in tail_copies:
                path_to_msd_csv = glob(os.path.join(path_to_processing_dir,
                                                      "data", "msd", f"{tail_copy}-{file_basename}"))[0]
                df_msd = pd.read_csv(path_to_msd_csv)
                if dt_ns is None:
                    dt_ns = df_msd["time_ns"].values[1] - df_msd["time_ns"].values[0]
                acorr_arrays.append(df_msd["msd"].values)

        time_ns, msd_avg = average_across_trajectories(acorr_arrays, dt_ns=dt_ns)
        df_acorr_avg = pd.DataFrame({"time_ns": time_ns, "msd": msd_avg})

        # Save the autocorrelation results to a CSV file
        os.makedirs(args.output_directory, exist_ok=True)
        df_acorr_avg.to_csv(os.path.join(args.output_directory, f"avg-{file_basename}"), index=False)
