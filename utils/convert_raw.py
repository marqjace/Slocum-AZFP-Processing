import os
import glob
import echopype as ep
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

if __name__ == "__main__":
    import sys
    convert_raw(sys.argv[1], sys.argv[2])