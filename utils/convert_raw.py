import os
import glob
import numpy as np
import echopype as ep
from scipy.interpolate import interp1d
import warnings

warnings.filterwarnings('ignore')

def calculate_azfp_sv_offset(pulse_length: float | int, frequency: float | int) -> float:
    """
    Created by Ian Black, PhD Student at Oregon State University, blackia@oregonstate.edu.

    Linearly interpolate an SV offset for a given pulse length and frequency.
    The predefined offsets used in this function can be found in Table 3 of the AZFP Operator's Manual.
    If the pulse length is smaller or larger than the pulse lengths provided in the manual, the offset is extrapolated using scipy.interpolate.interp1d. If the value is extrapolated, a UserWarning is raised.

    :param pulse_length: A pulse length value (in microseconds) from an AZFP config file for a given frequency channel.
    :param frequency: The AZFP frequency channel.
    :return: Either the known Sv offset if it is predefined, or an interpolated/extrapolated value for less common pulse lengths.
    """
    pulse_length = int(pulse_length) # Convert to an integer for consistency.
    frequency = int(frequency)

    high_frequencies = [120000,200000,455000,769000] # AZFP high frequency channels.
    sv_offset_hf = {  # Known high frequency channel correction values.
        150:1.4,
        200:1.4,
        250:1.3,
        300:1.2,
        500:0.9,
        700:0.6,
        900:0.3,
        1000:0.2
    }

    low_frequencies = [38000,67000,67500] # AZFP low frequency channels.
    sv_offset_lf = { # Known low frequency channel correction values.
        500:1.1,
        1000:0.7
    }

    if pulse_length < 150 or pulse_length > 1000:
        raise UserWarning("Are you sure your pulse length is less than 150 or greater than 1000 microseconds?")

    if frequency in high_frequencies:
        known_offsets = sv_offset_hf
    elif frequency in low_frequencies:
        if pulse_length < 300:
            raise UserWarning("Are you sure your pulse length is less than 300 microseconds?")
        known_offsets = sv_offset_lf

    if pulse_length in known_offsets.keys(): # If the pulse length already has a known offset, just return that without interpolating.
        known_offset = known_offsets[pulse_length]
        return known_offset
    else: # If the pulse length doesn't have a known offset, linearly interpolate.
        (x,y) = zip(*known_offsets.items())
        if min(x) <= pulse_length <= max(x):
            interpolated_offset = round(float(np.interp(pulse_length, x, y)),3)
        elif pulse_length < min(x) or pulse_length > max(x):
            func = interp1d(x,y, kind = 'linear', fill_value = 'extrapolate')
            interpolated_offset = round(float(func(pulse_length)),3)
        return interpolated_offset

pulse_length = 350
frequencies = [67000, 120000, 200000]

for frequency in frequencies:
    offset = calculate_azfp_sv_offset(pulse_length=pulse_length, frequency=frequency)
    ep.convert.parse_azfp.SV_OFFSET[frequency][pulse_length] = offset

def convert_raw(raw_data_directory, xml_file):
    """
    Convert raw AZFP echogram data to NetCDF using echopype.

    :param raw_data_directory: str, path to the raw echogram files
    :param xml_file: str, path to XML calibration file

    :return: ed_list (list of xarray datasets), proc_data_directory (path where NetCDF files are stored)
    """
    file_list = glob.glob(os.path.join(raw_data_directory, '*.01?'))
    file_list.sort()

    parent_path = os.path.dirname(os.path.normpath(raw_data_directory))
    proc_data_directory = os.path.join(parent_path, 'processed')
    os.makedirs(proc_data_directory, exist_ok=True)

    raw_directory = os.path.join(proc_data_directory, 'raw')
    os.makedirs(raw_directory, exist_ok=True)

    ed_list = []
    failed_files = []
    reused_files = []

    print(f"Processing {len(file_list)} AZFP binary files and converting them to NetCDF...")

    for i, raw_file in enumerate(file_list, 1):
        print(f"[{i}/{len(file_list)}] Processing file: {raw_file}")

        # Check if output already exists before parsing
        out_file = os.path.join(raw_directory, os.path.basename(raw_file) + '.raw.nc')
        
        if os.path.exists(out_file):
            print(f"Output file already exists, loading existing: {out_file}")
            try:
                ed = ep.open_converted(out_file)
                ed_list.append(ed)
                reused_files.append(raw_file)
            except Exception as e_load:
                print(f"Failed to load existing NetCDF {out_file}: {e_load}")
                failed_files.append(raw_file)
            continue

        try:
            ed = ep.open_raw(raw_file, sonar_model='AZFP', xml_path=xml_file)

            # Update platform metadata
            ed.platform.attrs['platform_name'] = 'Glider 592'
            ed.platform.attrs['platform_type'] = 'Sub-Surface Glider'
            ed.platform.attrs['platform_code_ICES'] = '27'

            # Replace NaN Sv_offset with zeros
            ed.vendor['Sv_offset'] = ed.vendor['DS'] * 0.0

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

    print(f"Finished processing {len(ed_list)} files total!")
    print(f"  - {len(ed_list) - len(reused_files)} newly converted")
    print(f"  - {len(reused_files)} reused from existing .raw.nc")
    
    if failed_files:
        print(f"\n{len(failed_files)} files failed during processing:")
        for f in failed_files:
            print(f"  {f}")

    return ed_list, proc_data_directory

if __name__ == "__main__":
    import sys
    convert_raw(sys.argv[1], sys.argv[2])