import scipy.io as sio
import xarray as xr
import numpy as np
from pathlib import Path
import numpy.ma as ma
import os

def convert_mat_to_netcdf(mat_file_path):
    """
    Converts a .mat file to a NetCDF file.

    :param mat_file_path: str, path to the .mat file to be converted
    :return: None, saves the converted NetCDF file in the same directory as the .mat file
    """
    matfile = sio.loadmat(mat_file_path)
    # for k, v in matfile.items():
        # if not k.startswith("__"):
        # #     print(f"{k}: shape={v.shape}, dtype={v.dtype}, first values={v.flat[:5]}")

        # # print(f"Loaded .mat file: {mat_file_path}")

    data_vars = {k: v for k, v in matfile.items() if not k.startswith("__")}
    xr_vars = {}

    for key, val in data_vars.items():
        val = np.squeeze(val)
        if np.issubdtype(val.dtype, np.floating):
            val = ma.masked_invalid(val)
        # Handle 2D arrays that look like date-time (N, 6)
        if val.ndim == 2 and val.shape[1] == 6 and np.issubdtype(val.dtype, np.integer):
            # Convert to datetime64[ns]
            dt = [np.datetime64(f"{row[0]:04d}-{row[1]:02d}-{row[2]:02d}T{row[3]:02d}:{row[4]:02d}:{row[5]:02d}")
                for row in val]
            xr_vars[key] = (["time"], np.array(dt))
        elif val.ndim == 1:
            xr_vars[key] = ([f"dim_{key}"], val)
        else:
            xr_vars[key] = (list([f"dim_{key}{i}" for i in range(val.ndim)]), val)
        
    ds_glider = xr.Dataset(xr_vars)
    print("Converted .mat file to xarray.Dataset")

    directory = Path(mat_file_path).parent
    filename = Path(mat_file_path).name
    name_only = filename.split('.')[0]
    nc_name = name_only + ".nc"

    ds_glider.to_netcdf(os.path.join(directory, nc_name), format='NETCDF4')
    # print(f"Saved NetCDF file: {os.path.join(directory, nc_name)}")

if __name__ == "__main__":
    import sys
    convert_mat_to_netcdf(sys.argv[1])