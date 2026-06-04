import glob
import os
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd
import cmocean
from dask.diagnostics import ProgressBar

import warnings
warnings.filterwarnings('ignore')

def merge_glider_azfp_dataset(proc_data_directory, transect_line=None):
    figures_directory = os.path.join(proc_data_directory, 'figures/')
    os.makedirs(figures_directory, exist_ok=True)

    # Define transect lines and corresponding time ranges
    if transect_line is None:
        transect_line = 'Complete Mission'
        start_time = None
        end_time = None
    elif transect_line == 1:
        transect_line = 'Line_1'
        start_time = pd.to_datetime('2023-05-24 18:20:00')
        end_time = pd.to_datetime('2023-05-26 07:38:00')
    elif transect_line == 2:
        transect_line = 'Line_2'
        start_time = pd.to_datetime('2023-05-26 07:38:00')
        end_time = pd.to_datetime('2023-05-28 23:22:00')
    elif transect_line == 3:
        transect_line = 'Line_3'
        start_time = pd.to_datetime('2023-05-28 23:22:00')
        end_time = pd.to_datetime('2023-05-30 17:53:00')
    elif transect_line == 4:
        transect_line = 'Line_4'
        start_time = pd.to_datetime('2023-05-30 17:53:00')
        end_time = pd.to_datetime('2023-06-02 10:27:00')
    elif transect_line == 5:
        transect_line = 'Line_5'
        start_time = pd.to_datetime('2023-06-02 10:27:00')
        end_time = pd.to_datetime('2023-06-04 02:51:00')
    elif transect_line == 6:
        transect_line = 'Line_6'
        start_time = pd.to_datetime('2023-06-04 02:51:00')
        end_time = pd.to_datetime('2023-06-05 16:15:00')

    file_list = glob.glob(os.path.join(proc_data_directory, '*_profiles.nc'))
    file_list.sort()

    print("Loading and merging processed AZFP data...")
    
    # Open lazily with optimized chunks
    ds_merged = xr.open_mfdataset(
        file_list,
        combine="nested",
        concat_dim="ping_time",
        coords="minimal",
        compat="override",
        parallel=True,
        engine="h5netcdf",
        chunks={"ping_time": 5000, "range_sample": -1}  # Keep range_sample together
    )
    
    # Filter BEFORE any computation
    if start_time is not None and end_time is not None:
        print(f"Filtering to time range: {start_time} to {end_time}")
        ds_merged = ds_merged.sel(ping_time=slice(start_time, end_time))
    
    # Save lazily (no .compute() needed - writes as it computes)
    output_file = os.path.join(proc_data_directory, f'sl592_{transect_line}.nc')
    print("Saving merged dataset...")
    with ProgressBar():
        ds_merged.to_netcdf(output_file, engine="h5netcdf")
    
    # NOW start plotting - compute only what's needed
    print("Preparing coordinate data for plotting...")
    
    # Compute small coordinate arrays once (these are lightweight)
    longitude = ds_merged["longitude"].values  # Small array, safe to compute
    depth = ds_merged["depth"].values
    bottom_depth = ds_merged["bottom_depth"].values
    channels = ds_merged["channel"].values
    
    # Plot echogram - compute per channel to save memory
    print("Plotting echogram...")
    fig, axes = plt.subplots(3, 1, figsize=(8, 12), dpi=150)
    
    for i, channel in enumerate(channels):
        print(f"  Loading channel {channel}...")
        # Compute ONLY this channel's data
        sv = ds_merged["Sv"].sel(channel=channel).values  # Auto-computes when .values called
        
        pcm = axes[i].pcolormesh(
            longitude, depth, sv,
            shading='auto', cmap='jet',
            vmin=-90, vmax=-60
        )
        axes[i].plot(longitude, bottom_depth, 'k-', linewidth=3)
        axes[i].invert_yaxis()
        axes[i].set_ylim(depth.max(), 0)
        axes[i].set_ylabel("Depth (m)")
        fig.colorbar(pcm, ax=axes[i], 
                    label="Volume Backscattering Strength\n(Sv re 1 $m^{-1}$) [dB]")
        
        # sv goes out of scope here, memory can be freed
    
    axes[0].set_title('67 kHz')
    axes[1].set_title('125 kHz')
    axes[2].set_title('200 kHz')
    plt.xlabel("Longitude ($\degree$E)")
    plt.suptitle(f"{transect_line}")
    plt.tight_layout()
    plt.savefig(
        os.path.join(figures_directory, f'{transect_line}_echogram.png'),
        dpi=150, bbox_inches='tight'
    )
    plt.close(fig)
    
    # Plot variables efficiently with config dictionary
    print("Plotting environmental variables...")
    var_configs = {
        'temperature': {
            'cmap': cmocean.cm.thermal,
            'label': 'Temperature (°C)',
            'vmin': 5, 'vmax': 14
        },
        'salinity': {
            'cmap': cmocean.cm.haline,
            'label': 'Salinity (PSU)',
            'vmin': 30, 'vmax': 35
        },
        'backscatter': {
            'cmap': cmocean.cm.algae,
            'label': 'Backscatter ($m^{-1}$)',
            'vmin': 0, 'vmax': 0.02
        },
        'chlorophyll': {
            'cmap': cmocean.cm.algae,
            'label': 'Chlorophyll (µg/L)',
            'vmin': 0, 'vmax': 15
        },
        'cdom': {
            'cmap': cmocean.cm.algae,
            'label': 'CDOM (ppb/L)',
            'vmin': 0, 'vmax': 6
        },
        'oxygen': {
            'cmap': cmocean.cm.oxy,
            'label': 'Oxygen (µmol/kg)',
            'vmin': None, 'vmax': None
        }
    }
    
    for var, config in var_configs.items():
        print(f"  Plotting {var}...")
        
        # Compute ONLY this variable's data
        var_data = ds_merged[var].values  # Auto-computes when .values called
        
        fig, ax = plt.subplots(figsize=(12, 6))
        pcm = ax.pcolormesh(
            longitude, depth, var_data,
            shading='auto',
            cmap=config['cmap'],
            vmin=config['vmin'],
            vmax=config['vmax']
        )
        ax.plot(longitude, bottom_depth, 'k-', linewidth=3)
        ax.invert_yaxis()
        ax.set_ylim(depth.max(), 0)
        ax.set_ylabel("Depth (m)")
        ax.set_xlabel("Longitude ($\degree$E)")
        ax.set_title(f"{transect_line}")
        fig.colorbar(pcm, ax=ax, label=config['label'])
        plt.tight_layout()
        plt.savefig(
            os.path.join(figures_directory, f'{var}_{transect_line}.png'),
            dpi=150, bbox_inches='tight'
        )
        plt.close(fig)
            
    ds_merged.close()
    print("Done!")


proc_data_directory = r"C:\Users\marqjace\data\azfp\processed\proc"

merge_glider_azfp_dataset(proc_data_directory)