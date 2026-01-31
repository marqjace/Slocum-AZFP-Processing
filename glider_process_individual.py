# glider_process_individual.py

# Created by Jace Marquardt
# Last updated: 2025-08-27

import echopype as ep
import os
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd
from convert_mat_to_netcdf import convert_mat_to_netcdf
from pathlib import Path
import matplotlib.dates as mdates
from convert_raw import convert_raw

import warnings
warnings.filterwarnings('ignore')

def glider_process_individual(data_directory, xml_file, glider_data):
    """
    This function takes the converted and calibrated echogram data from the glider AZFP, processes it, and generates echograms.
    It also saves the echograms as PNG files in the specified figures directory within '/processed'.

    :param data_directory: str, path to the directory containing the raw echogram data files
    :param proc_data_directory: str, path to the directory where the converted netCDF files will be saved
    :param xml_file: str, path to the XML file containing the sonar configuration
    :param glider_data: matlab file
    
    :return: Saved echograms as PNG files in the figures directory
    """
    ############################ Glider Data Processing ############################

    # Convert the .mat glider data file to a netCDF file
    convert_mat_to_netcdf(glider_data)
    print('Converting the glider data from MatLab to NetCDF....')

    directory = Path(glider_data).parent
    filename = Path(glider_data).name
    name_only = filename.split('.')[0]
    nc_name = name_only + ".nc"
    
    # Load the glider data as an xarray dataset
    with xr.open_dataset(os.path.join(directory, nc_name)) as ds_glider:
        print(f"Loaded glider data from: {os.path.join(directory, nc_name)}")

        # Extract temperature, salinity, and depth from the glider dataset
        temperature = ds_glider['Temp'].values
        salinity = ds_glider['Salt'].values
        depth = ds_glider['Depth'].values
        bottom_depth = ds_glider['Bottom_depth'].values
        glider_time = ds_glider['Timeinsec'].values
        back_scatter = ds_glider['BackScatter'].values
        cdom = ds_glider['Cdom'].values
        chlor = ds_glider['Chlor'].values
        oxygen = ds_glider['Oxygen'].values
        lon = ds_glider['Lon'].values
        lat = ds_glider['Lat'].values

    # Mask the data where there are NaNs
    valid = (~np.isnan(depth)) & (~np.isnan(glider_time))
    glider_time_valid = glider_time[valid]
    temperature_valid = temperature[valid]
    salinity_valid = salinity[valid]
    depth_valid = depth[valid]
    bottom_depth_valid = bottom_depth[valid]
    back_scatter_valid = back_scatter[valid]
    cdom_valid = cdom[valid]
    chlor_valid = chlor[valid]
    oxygen_valid = oxygen[valid]
    lon_valid = lon[valid]
    lat_valid = lat[valid]

    # Sort by time
    sort_idx = np.argsort(glider_time_valid)
    glider_time_valid = glider_time_valid[sort_idx]
    depth_valid = depth_valid[sort_idx]
    bottom_depth_valid = bottom_depth_valid[sort_idx]
    temperature_valid = temperature_valid[sort_idx]
    salinity_valid = salinity_valid[sort_idx]
    back_scatter_valid = back_scatter_valid[sort_idx]
    cdom_valid = cdom_valid[sort_idx]
    chlor_valid = chlor_valid[sort_idx]
    oxygen_valid = oxygen_valid[sort_idx]
    lon_valid = lon_valid[sort_idx]
    lat_valid = lat_valid[sort_idx]

    ############################ AZFP Data Processing ############################

    # Load the AZFP data
    ed_list, proc_data_directory = convert_raw(data_directory, xml_file)
    print('Processing raw AZFP files....')
    
    processed_directory = os.path.join(proc_data_directory, 'proc')
    if not os.path.isdir(processed_directory):
        os.mkdir(processed_directory)
    
    print(f"Processing {len(ed_list)} raw AZFP NetCDF files...")

    for i, ds in enumerate(ed_list):
        print(f'Processing file: ({i+1}/{len(ed_list)})')
        file = os.path.basename(ds.provenance.source_filenames.values[0])

        # Interpolate glider data
        ping_time = ds.environment.time1.values.astype(float)
        t = np.interp(ping_time, glider_time_valid, temperature_valid)
        s = np.interp(ping_time, glider_time_valid, salinity_valid)
        d = np.interp(ping_time, glider_time_valid, depth_valid)

        # print(np.nanmean(t), np.nanmean(s), np.nanmean(d))

        # Set the environmental parameters
        env_params = {
            'pressure': np.nanmean(d),
            'temperature': np.nanmean(t),
            'salinity': np.nanmean(s)
        }

        ds_sv = ep.calibrate.compute_Sv(ds, env_params=env_params) # Compute Sv
        # ds_sv = ds_sv.isel(range_sample=slice(35, None)) # Remove near-surface data (first 35 samples)
        ds_sv_clean = ep.clean.remove_background_noise(ds_sv, ping_num=5, range_sample_num=30) # Remove background noise
        ds_sv.close()

        ############################### AZFP Time Correction ###############################
        # print("Correcting AZFP ping times...")

        # Interpolate AZFP time to match glider time and convert to julian date
        ds_sv_clean['ping_time'] = pd.to_datetime(ds_sv_clean['ping_time'].values, unit='s') + pd.to_timedelta(505, unit='s')
        ping_time = ds_sv_clean['ping_time'].values

        ############################### Mask Glider Time to AZFP Time Range ###############################
        # Convert glider time to julian date
        glider_dt = pd.to_datetime(glider_time_valid, unit='s').values

        # Determine the AZFP time range
        start_julday = ping_time.min()
        end_julday = ping_time.max()

        # Mask the glider time to AZFP time range
        mask = ~np.isnan(depth_valid) & ~np.isnan(glider_dt) & (glider_dt >= start_julday) & (glider_dt <= end_julday)
        glider_dt_masked = glider_dt[mask]
        depth_masked = depth_valid[mask]
        bottom_depth_masked = bottom_depth_valid[mask]
        temperature_masked = temperature_valid[mask]
        salinity_masked = salinity_valid[mask]
        back_scatter_masked = back_scatter_valid[mask]
        cdom_masked = cdom_valid[mask]
        chlor_masked = chlor_valid[mask]
        oxygen_masked = oxygen_valid[mask]
        lon_masked = lon_valid[mask]
        lat_masked = lat_valid[mask]

        lon_mean = np.nanmean(lon_masked)
        lat_mean = np.nanmean(lat_masked)
        bottom_depth_mean = np.nanmean(bottom_depth_masked)

        ############################### Glider Depth onto AZFP Ping Times ###############################
        # print("Interpolating glider depth onto AZFP ping times...")

        # Interpolate glider depth only onto those times
        glider_depth_interp = np.interp(
            mdates.date2num(ping_time),             
            mdates.date2num(glider_dt_masked),      
            depth_masked                       
        )

        # Convert glider depth to a DataArray with ping_time coordinates
        glider_depth_da = xr.DataArray(
            glider_depth_interp,
            dims=["ping_time"],
            coords={"ping_time": ds_sv_clean["ping_time"]}
        )

        # Expand it to 2D to match echo_range
        glider_depth_2d = glider_depth_da.broadcast_like(ds_sv_clean["echo_range"])

        # Add to echo_range
        echo_range_depth = ds_sv_clean["echo_range"] + glider_depth_2d

        # Reassign echo_range to echo_range_depth values
        ds_sv_clean = ds_sv_clean.assign_coords(
            echo_range=echo_range_depth
        )

        mask_below_bottom = ds_sv_clean["echo_range"] > bottom_depth_mean
        ds_sv_clean = ds_sv_clean.where(~mask_below_bottom)

        # # ############################### Masking step ###############################
        # # # Perform bottom detection and masking before saving
        # # channels = ds_sv_clean["channel"].values
        # threshold = -45.0

        ################################ Masking + Depth-binned profiles in one loop ################################
        sv_profiles = []
        glider_profiles = []

        for channel in ds_sv_clean["channel"].values:
            sv = ds_sv_clean["Sv_corrected"].sel(channel=channel)

            # ---------------- Mask out climb ----------------
            mean_depth = ds_sv_clean['echo_range'].mean(
                dim=[d for d in ds_sv_clean['echo_range'].dims if d != "ping_time"]
            )
            depth_derivative = mean_depth.diff("ping_time")
            depth_derivative = depth_derivative.reindex(
                {"ping_time": ds_sv_clean["ping_time"]}, method="nearest"
            )
            climb_mask = depth_derivative < 0
            sv = sv.where(~climb_mask)

            # # ---------------- Uncomment to Plot Individual Dive Echograms ----------------
            # echo_range = ds_sv_clean["echo_range"].sel(channel=channel)
            # C = sv.values.T  # shape (M, N)
            # Y = echo_range.transpose("range_sample", "ping_time").values  # (M, N)
            # X = np.broadcast_to(sv['ping_time'].values.reshape(1, -1), C.shape)  # (M, N)

            # fig, ax = plt.subplots(figsize=(12, 6))

            # pcolormesh = ax.pcolormesh(
            #     X,
            #     Y,
            #     C,
            #     shading='auto',
            #     vmin=-100,
            #     vmax=-60,
            #     cmap='jet',
            # )

            # ax.xaxis.set_major_locator(mdates.AutoDateLocator())
            # ax.xaxis.set_major_formatter(mdates.DateFormatter('%d %H:%M'))

            # # Rotate labels for readability
            # fig.autofmt_xdate(rotation=30, ha='right')

            # # ax.invert_yaxis()
            # ax.set_ylim(Y.max(), 0)
            # ax.set_title(f"{channel} Sv")
            # ax.set_ylabel("Depth (m)")
            # ax.set_xlabel("Ping Time (DD HH:MM)")
            # fig.colorbar(pcolormesh, ax=ax, label="Volume backscattering strength (Sv re 1 m-1) [dB]")

            # plt.savefig(os.path.join(processed_directory, f'{file}_{channel}_echogram.png'), dpi=300, bbox_inches='tight')
            # plt.close(fig)

            # ---------------- Depth-binning Sv ----------------
            depth = ds_sv_clean["echo_range"].sel(channel=channel)
            depth_flat = depth.values.ravel()
            sv_flat = sv.values.ravel()

            mask = np.isfinite(depth_flat) & np.isfinite(sv_flat)
            depth_flat = depth_flat[mask]
            sv_flat = sv_flat[mask]

            if depth_flat.size == 0 or sv_flat.size == 0:
                print(f"Skipping {file} {channel}: no valid data")
                continue

            depth_bins = np.arange(0, np.nanmax(depth_flat) + 1, 1)
            bin_indices = np.digitize(depth_flat, depth_bins)

            sv_binned = [
                np.nanmean(sv_flat[bin_indices == i]) if np.any(bin_indices == i) else np.nan
                for i in range(1, len(depth_bins))
            ]

            mean_ping_time = pd.to_datetime(ds_sv_clean["ping_time"].values).mean()

            sv_profiles.append(
                xr.DataArray(
                    np.array(sv_binned)[:, np.newaxis],
                    dims=["depth", "ping_time"],
                    coords={
                        "depth": depth_bins[:-1],
                        "ping_time": [mean_ping_time],
                        "longitude": (["ping_time"], [lon_mean]),   # <-- associate with ping_time
                        "latitude": (["ping_time"], [lat_mean]),    # <-- associate with ping_time
                        "bottom_depth": (["ping_time"], [bottom_depth_mean]),  # <-- associate with ping_time
                    },
                    name="Sv"
                )
            )

        # ---------------- Concat Sv across channel ----------------
        if len(sv_profiles) > 0:
            sv_profiles_ds = xr.concat(sv_profiles, dim="channel")
            sv_profiles_ds = sv_profiles_ds.assign_coords(channel=ds_sv_clean["channel"].values)
        else:
            sv_profiles_ds = None

        # ---------------- Depth-binned Glider Variables (only once, not per channel) ----------------
        glider_vars = {
            "temperature": temperature_masked,
            "salinity": salinity_masked,
            "backscatter": back_scatter_masked,
            "cdom": cdom_masked,
            "chlorophyll": chlor_masked,
            "oxygen": oxygen_masked,
        }

        for var_name, var_values in glider_vars.items():
            var_interp = np.interp(
                mdates.date2num(ping_time),
                mdates.date2num(glider_dt_masked),
                var_values
            )

            var_binned = [
                np.nanmean(var_interp[(glider_depth_interp >= depth_bins[i-1]) &
                                    (glider_depth_interp < depth_bins[i])])
                if np.any((glider_depth_interp >= depth_bins[i-1]) &
                        (glider_depth_interp < depth_bins[i])) else np.nan
                for i in range(1, len(depth_bins))
            ]

            glider_profiles.append(
                xr.DataArray(
                    np.array(var_binned)[:, np.newaxis],
                    dims=["depth", "ping_time"],
                    coords={
                        "depth": depth_bins[:-1],
                        "ping_time": [mean_ping_time],
                        "longitude": (["ping_time"], [lon_mean]),   # <-- associate with ping_time
                        "latitude": (["ping_time"], [lat_mean]),    # <-- associate with ping_time
                        "bottom_depth": (["ping_time"], [bottom_depth_mean]),  # <-- associate with ping_time
                    },
                    name=var_name
                )
            )

        glider_profiles_ds = xr.merge(glider_profiles) if glider_profiles else None

        # ---------------- Merge Sv + glider vars ----------------
        if sv_profiles_ds is not None and glider_profiles_ds is not None:
            profiles_ds = xr.merge([sv_profiles_ds, glider_profiles_ds])
        elif sv_profiles_ds is not None:
            profiles_ds = sv_profiles_ds
        elif glider_profiles_ds is not None:
            profiles_ds = glider_profiles_ds
        else:
            print(f"Skipping {file}: no valid profiles")
            continue

        profiles_ds.to_netcdf(os.path.join(processed_directory, f"{file}_profiles.nc"))
        profiles_ds.close()


# glider_data = r"C:\Users\marqjace\azfp\2023_deployment\processing\WA_202305241820-deployment_osu592_pass3.mat"
# raw_data_directory = r"C:\Users\marqjace\azfp\2023_deployment\processing\to_process"
xml_file = 'tweaked.xml'
# xml_file = 'tweaked2.xml'

glider_data = r"C:\Users\marqjace\azfp\2023_deployment\processing\WA_202305241820-deployment_osu592_pass3.mat"
raw_data_directory = r"C:\Users\marqjace\azfp\2023_deployment\processing\to_process"

glider_process_individual(raw_data_directory, xml_file, glider_data)