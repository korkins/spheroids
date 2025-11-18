import os
import numpy as np

folder1 = "./github_bins_bmrk"
folder2 = "./KERNEL_n22_fix"

files = sorted(os.listdir(folder1))

for fname in files:
    f1 = os.path.join(folder1, fname)
    f2 = os.path.join(folder2, fname)

    if not os.path.isfile(f2):
        print(f"{fname}: NOT FOUND in folder2")
        continue

    a = np.fromfile(f1, dtype=np.float64)
    b = np.fromfile(f2, dtype=np.float64)

    if a.size != b.size:
        print(f"{fname}: size mismatch ({a.size} vs {b.size})")
        continue

    maxdiff = np.max(np.abs(a - b))
    print(f"{fname}: max abs diff = {maxdiff:.1e}")