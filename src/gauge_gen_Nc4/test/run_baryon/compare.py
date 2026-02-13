#!/usr/bin/env python3

import h5py
import numpy as np
import sys


# to read h5data with re/im
dt=np.dtype([('re',np.double),('im',np.double)])

def compare_datasets(dset1, dset2, path):
    # print(dset1,dset2)
    if dset1.shape != dset2.shape:
        print(f"[{path}] Shape mismatch: {dset1.shape} vs {dset2.shape}")
        return False

    diff = np.abs(dset1[...] - dset2[...])
    # print(dset1[...] - dset2[...]) 
    if not np.allclose(dset1[...], dset2[...], rtol=1e-10, atol=1e-12):
        print(f"[{path}] Values differ (max diff = {diff.max()})")
        return False
    return True

def compare_groups(g1, g2, path="/"):
    ok = True
    for name in g1:
        print(name)
        new_path = f"{path}{name}/"
        if name not in g2:
            print(f"{new_path} missing in second file")
            ok = False
            continue

        tmp=np.array(g1[name]['data'][:],dtype=dt)
        item1=tmp['re']+1j*tmp['im']
        tmp=np.array(g2[name]['data'][:],dtype=dt)
        item2=tmp['re']+1j*tmp['im']
        # f=g1[name]['data'][:]
        # print(f)
        # print(g2)

        if not compare_datasets(item1, item2, new_path):
            ok = False
            
    return ok

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python compare_h5.py file1.h5 file2.h5")
        sys.exit(1)

    file1, file2 = sys.argv[1], sys.argv[2]
    with h5py.File(file1, 'r') as f1, h5py.File(file2, 'r') as f2:
        result = compare_groups(f1, f2)
        if result:
            print("✅ Files match numerically.")
        else:
            print("❌ Files differ.")
