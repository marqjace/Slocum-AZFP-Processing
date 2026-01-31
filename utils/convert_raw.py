# # convert_raw.py
# # Converting the raw data from the glider AZFP to netCDF format using echopype

# # Created by Jace Marquardt
# # Last updated: 2025-06-12

# import echopype as ep
# import glob
# import os

# import warnings
# warnings.filterwarnings('ignore')

# def convert_raw(data_directory, xml_file):
#     """
#     This function takes the raw echogram data from the glider AZFP, calibrates it using the XML file, and converts it to netCDF format.

#     :param data_directory: str, path to the directory containing the raw echogram data files
#     :param xml_file: str, path to the XML file containing the sonar configuration

#     :return: Converted and calibrated bioacoustic sonar data in a xarray dataset object
#     """

#     file_list = glob.glob(os.path.join(data_directory, '*.01?'))
#     file_list.sort()

#     parent_path = os.path.dirname(data_directory)

#     proc_data_directory = os.path.join(parent_path, 'processed')
#     if not os.path.isdir(proc_data_directory):
#         os.mkdir(proc_data_directory)

#     raw_directory = os.path.join(proc_data_directory, 'raw')
#     if not os.path.isdir(raw_directory):
#         os.mkdir(raw_directory)

#     # Convert the list of .01X files using echopype and save the output as netCDF files
#     ed_list = []

#     print(f"Processing {len(file_list)} AZFP binary files and converting them to NetCDF...")

#     for i, raw_file in enumerate(file_list):
#         print(f'Processing file: {raw_file} ({i+1}/{len(file_list)})')

#         try:
#             # Load the raw data into an xarray dataset
#             ed = ep.open_raw(raw_file, sonar_model='AZFP', xml_path=xml_file)
#         except Exception as e:
#             print(f"Skipping file due to error: {raw_file}")
#             print(f"   Error: {e}")
#             continue  # Skip to next file

#         # Update platform metadata
#         ed.platform.attrs['platform_name'] = 'Glider 592'
#         ed.platform.attrs['platform_type'] = 'Sub-Surface Glider'
#         ed.platform.attrs['platform_code_ICES'] = '27'

#         # Change Sv_offset values from NaN to 0's
#         ed.vendor['Sv_offset'] = ed.vendor['DS'] * 0.0

#         # Save the converted data to netCDF
#         ed.to_netcdf(os.path.join(raw_directory, os.path.split(raw_file)[1] + '.raw.nc'))
#         ed_list.append(ed)

#     print(f"Finished processing {len(ed_list)} AZFP binary files!")

#     return ed_list, proc_data_directory

# # raw_data_directory = r"C:\Users\marqjace\azfp\2024_deployment\processing\to_process"
# # xml_file = 'tweaked_2024.xml'

# # ed_list, proc_data_directory = convert_raw(raw_data_directory, xml_file)

# # print(ed_list[0])





# convert_raw.py
# Converting the raw data from the glider AZFP to netCDF format using echopype

import echopype as ep
import glob
import os
import warnings

warnings.filterwarnings('ignore')

def convert_raw(data_directory, xml_file):
    """
    Convert raw AZFP echogram data to NetCDF using echopype.

    :param data_directory: str, path to the raw echogram files
    :param xml_file: str, path to XML calibration file

    :return: ed_list (list of xarray datasets), proc_data_directory (path where NetCDF files are stored)
    """
    file_list = glob.glob(os.path.join(data_directory, '*.01?'))
    file_list.sort()

    parent_path = os.path.dirname(data_directory)
    proc_data_directory = os.path.join(parent_path, 'processed')
    os.makedirs(proc_data_directory, exist_ok=True)

    raw_directory = os.path.join(proc_data_directory, 'raw')
    os.makedirs(raw_directory, exist_ok=True)

    ed_list = []
    failed_files = []

    print(f"Processing {len(file_list)} AZFP binary files and converting them to NetCDF...")

    for i, raw_file in enumerate(file_list, 1):
        print(f"[{i}/{len(file_list)}] Processing file: {raw_file}")

        try:
            ed = ep.open_raw(raw_file, sonar_model='AZFP', xml_path=xml_file)

            # Update platform metadata
            ed.platform.attrs['platform_name'] = 'Glider 592'
            ed.platform.attrs['platform_type'] = 'Sub-Surface Glider'
            ed.platform.attrs['platform_code_ICES'] = '27'

            # Replace NaN Sv_offset with zeros
            ed.vendor['Sv_offset'] = ed.vendor['DS'] * 0.0

            # Save to NetCDF with safe format
            out_file = os.path.join(raw_directory, os.path.basename(raw_file) + '.raw.nc')

            # Avoid overwriting old files
            if os.path.exists(out_file):
                print(f"Output file already exists, skipping: {out_file}")
                failed_files.append(raw_file)
                continue

            try:
                ed.to_netcdf(out_file, format='NETCDF4')
            except Exception as e_nc:
                print(f"Failed to save NetCDF for {raw_file}, skipping.")
                print(f"Error: {e_nc}")
                failed_files.append(raw_file)
                continue

            ed_list.append(ed)

        except Exception as e:
            print(f"Failed to process raw file {raw_file}, skipping.")
            print(f"Error: {e}")
            failed_files.append(raw_file)
            continue

    print(f"Finished processing {len(ed_list)} files!")
    if failed_files:
        print(f"{len(failed_files)} files failed during processing:")
        for f in failed_files:
            print(f"{f}")

    return ed_list, proc_data_directory

