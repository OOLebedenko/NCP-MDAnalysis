import os
import argparse
from glob import glob

import pandas as pd
from tqdm import tqdm

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Calculation N-H autocorrelation function averaged over trajectories of the same length')
    parser.add_argument('--path-to-processing-dirs', required=True)
    parser.add_argument('--output-directory', default=".")
    args = parser.parse_args()

    path_to_processing_dirs = [path_to_processing_dir.lstrip().strip()
                               for path_to_processing_dir in args.path_to_processing_dirs.split(",")]
    vector_basenames = [os.path.basename(acorr_csv).split("-", 1)[-1] for acorr_csv in
                        glob(os.path.join(path_to_processing_dirs[0], "H4_1", "data", "acorr", "*.csv"))]
    vector_basenames.sort()

    tail_copies = ["H4_1", "H4_2"]
    time_ns = None
    for vector_basename in tqdm(vector_basenames, desc="average autocorrelation function"):
        df_acorr_avg = None
        for path_to_processing_dir in path_to_processing_dirs:
            for tail_copy in tail_copies:
                path_to_acorr_csv = glob(os.path.join(path_to_processing_dir, tail_copy,
                                                      "data", "acorr", f"*{vector_basename}"))[0]
                df_acorr = pd.read_csv(path_to_acorr_csv)

                if df_acorr_avg is None:
                    df_acorr_avg = df_acorr
                else:
                    df_acorr_avg += df_acorr

        df_acorr_avg = df_acorr_avg / len(path_to_processing_dirs) / len(tail_copies)

        # Save the autocorrelation results to a CSV file
        os.makedirs(args.output_directory, exist_ok=True)
        df_acorr_avg.to_csv(os.path.join(args.output_directory, f"avg-{vector_basename}"), index=False)
