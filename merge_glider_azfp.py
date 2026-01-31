import glob
import os
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd
import cmocean


import warnings
warnings.filterwarnings('ignore')

def merge_glider_azfp_dataset(proc_data_directory, transect_line=None):
        
        # Make figures directory
        figures_directory = os.path.join(proc_data_directory, 'figures/')
        if not os.path.isdir(figures_directory):
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

        # Concatenate along 'ping_time' without fully loading into memory
        print("Loading and merging processed AZFP data...")

        # ds_merged = xr.open_mfdataset(
        #     file_list,
        #     combine="by_coords",
        #     join="outer",
        #     compat="override",
        #     coords="minimal",
        #     parallel=False,
        #     chunks={}
        # )

        ds_merged = xr.open_mfdataset(
            file_list,
            combine="nested",
            concat_dim="ping_time",
            coords="minimal",
            compat="override",
            parallel=True,
            engine="h5netcdf",
            chunks={"ping_time": 1000}  # tune based on dataset size
        )

        print("Loading data into memory...")
        ds_merged = ds_merged.compute()  # trigger computation now

        # print("Saving merged dataset...")
        # ds_merged.to_netcdf(
        #     os.path.join(proc_data_directory, 'sl592_complete_azfp.nc'),
        #     engine="netcdf4"  # or "h5netcdf" if you prefer
        # )

        print("Saving merged dataset...")
        # ds_merged.to_netcdf(os.path.join(proc_data_directory, 'sl592_complete_azfp.nc'))
        ds_merged.to_netcdf(
            os.path.join(proc_data_directory, 'sl592_complete_azfp_01_30_2026.nc'),
            engine="h5netcdf"
        )

        # # Mask specified time range
        # ds_merged = ds_merged.sel(ping_time=slice(start_time, end_time))

        # channels = ds_merged["channel"].values

        # # Create figure
        # print("Plotting merged echogram...")
        # fig, axes = plt.subplots(3,1, figsize=(8, 12), dpi=300)

        # for i, channel in enumerate(channels):
        #     # print(f"Plotting for channel: {channel}")
        #     sv = ds_merged["Sv"].sel(channel=channel)
        #     longitude = ds_merged["longitude"].values
        #     depth = ds_merged["depth"].values

        #     pcm = axes[i].pcolormesh(
        #         longitude,
        #         depth,
        #         sv.values,
        #         shading='auto',
        #         cmap='jet',
        #         vmin=-90,
        #         vmax=-60
        #     )

        #     bottom = axes[i].plot(longitude, ds_merged["bottom_depth"].values, 'k-', linewidth=5)

        #     axes[i].invert_yaxis()
        #     axes[i].set_ylim(depth.max(), 0)
        #     axes[i].set_ylabel("Depth (m)")
        #     axes[0].set_title('67 kHz')
        #     axes[1].set_title('125 kHz')
        #     axes[2].set_title('200 kHz')
            
        #     fig.colorbar(pcm, ax=axes[i], label=f"Volume Backscattering Strength\n(Sv re 1 $m^-1$) [dB]")

        # plt.xlabel("Longitude ($\degree$E)")
        # plt.suptitle(f"{transect_line}, Time (UTC): {str(start_time)} - {str(end_time)}")
        # plt.tight_layout()
        # plt.savefig(os.path.join(figures_directory, f'{transect_line}_echogram.png'), dpi=300)
        # plt.close(fig)

        # # ------------------------ Plot Glider Variables ------------------------
        # variables = ['temperature', 'salinity', 'backscatter', 'cdom', 'chlorophyll', 'oxygen']

        # for var in variables:
        #     print(f"Plotting {var}")
        #     # Create figure
        #     fig, ax = plt.subplots(figsize=(12, 6))

        #     # Convert ping_time to matplotlib dates
        #     longitude = ds_merged["longitude"].values
        #     depth = ds_merged["depth"].values

        #     if var == 'temperature':
        #         cmap = cmocean.cm.thermal
        #         var_label = 'Temperature (°C)'
        #         vmin = 5
        #         vmax = 14

        #     elif var == 'salinity':
        #         cmap = cmocean.cm.haline
        #         var_label = 'Salinity (PSU)'
        #         vmin = 30
        #         vmax = 35

        #     elif var == 'backscatter':
        #         cmap = cmocean.cm.algae
        #         var_label = f'Backscatter ($m^-1$)'
        #         vmin = 0
        #         vmax = 0.02

        #     elif var == 'chlorophyll':
        #         cmap = cmocean.cm.algae
        #         var_label = 'Chlorophyll (µg/L)'
        #         vmin = 0
        #         vmax = 15

        #     elif var == 'cdom':
        #         cmap = cmocean.cm.algae
        #         var_label = 'CDOM (ppb/L)'
        #         vmin = 0
        #         vmax = 6

        #     elif var == 'oxygen':
        #         cmap = cmocean.cm.oxy
        #         var_label = 'Oxygen (µmol/kg)'
        #         vmin = None
        #         vmax = None

        #     else:
        #         cmap = 'viridis'
        #         var_label = var
        #         vmin = None
        #         vmax = None

        #     # Plot pcolormesh
        #     pcm = ax.pcolormesh(
        #         longitude,
        #         depth,
        #         ds_merged[var].values,
        #         shading='auto',
        #         cmap=cmap,
        #         vmin=vmin,
        #         vmax=vmax
        #     )

        #     bottom = ax.plot(longitude, ds_merged["bottom_depth"].values, 'k-', linewidth=5)

        #     # Format axes
        #     ax.invert_yaxis()
        #     ax.set_ylim(depth.max(), 0)
        #     ax.set_ylabel("Depth (m)")
        #     ax.set_xlabel(f"Longitude ($\degree$E)")
        #     ax.set_title(f"{transect_line}, Time (UTC): {str(start_time)} - {str(end_time)}")
        #     fig.colorbar(pcm, ax=ax, label=f"{var_label}")
        #     plt.tight_layout()
        #     plt.savefig(os.path.join(figures_directory, f'{var}_{transect_line}.png'), dpi=300, bbox_inches='tight')
        #     plt.close(fig)

        ds_merged.close()
        print("Done!")


proc_data_directory = r"C:\Users\marqjace\azfp\2023_deployment\processing\processed\proc"

merge_glider_azfp_dataset(proc_data_directory)