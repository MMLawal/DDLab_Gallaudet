import numpy as np
import pandas as pd
import glob

results = []

for file in glob.glob("GIST_grid_*.npy"):

    grid = np.load(file)

    isoform = file.split("_")[2]
    pocket = file.split("_")[3]

    values = grid[~np.isnan(grid)]

    mean_dG = np.mean(values)

    unstable_fraction = np.sum(values > 0) / len(values)

    entropy_proxy = np.var(values)

    results.append({
        "isoform": isoform,
        "pocket": pocket,
        "mean_dG": mean_dG,
        "fraction_unstable": unstable_fraction,
        "hydration_entropy": entropy_proxy
    })

df = pd.DataFrame(results)

pocket_summary = df.groupby(["isoform","pocket"]).mean()

print(pocket_summary)
