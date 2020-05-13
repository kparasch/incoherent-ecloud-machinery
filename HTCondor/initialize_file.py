import h5py
import sys

out_file = sys.argv[1]
with h5py.File(out_file, 'w-') as fid:
    grp = fid.create_group('checkpoint')
    grp['checkpoint'] = 0
