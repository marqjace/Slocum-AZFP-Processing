# Slocum-AZFP-Processing

## glider_process_individual.py
To run: `uv run glider_process_individual.py`
- Edit lines 30-32 for correct filepaths
- Uncomment lines 224-256 to create Sv figures for each file

Utilizes **convert_raw.py** and **convert_mat_to_netcdf.py** under `\utils`. 
- **convert_raw.py** converts raw AZFP echogram data to NetCDF using echopype.
- **convert_mat_to_netcdf.py** converts a .mat file to a NetCDF file.

## merge_glider_azfp.py
To run: `uv run merge_glider_azfp.py`
- Edit lines 17-44 for correct transect line dates