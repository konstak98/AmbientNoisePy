import os
from obspy import read, UTCDateTime
from obspy.io.sac import SACTrace
from obspy.core import AttribDict
from scipy.signal import butter, filtfilt
from obspy.signal.invsim import cosine_taper
import numpy as np
from datetime import datetime

# Function to whiten the input trace
def whiten_trace(trace, freqmin, freqmax):
    # Implementation details omitted for brevity...
    pass

# Function to get filenames based on file type and extension
def get_filenames(file_type, extension, working_directory):
    if file_type == 2:
        extension_with_dot = f".{extension}"
        all_files_paths = [os.path.join(root, file) for root, dirs, files in os.walk(working_directory) for file in files if file.endswith(extension_with_dot)]
    else:
        extension_lower = f".r{extension[1]}{extension[2]}"
        extension_upper = f".R{extension[1]}{extension[2]}"
        all_files_paths = [os.path.join(root, file) for root, dirs, files in os.walk(working_directory) for file in files if file.endswith(extension_lower) or file.endswith(extension_upper)]
    return all_files_paths

# Function to retrieve station coordinates
def get_station_coordinates(delimiter, coordinates_file_path, station_id):
    with open(coordinates_file_path) as file:
        for line in file:
            parts = line.split(delimiter)
            current_station_id = parts[0]
            latitude = parts[1]
            longitude = parts[2]
            if current_station_id == station_id:
                return latitude, longitude

    # If station ID is not found, return None
    return None, None

def process_mseed_files(input_mseed_directory, input_sac_directory, full_coord_file_path, channels):
    file_paths = get_filenames(2, 'mseed', input_mseed_directory)

    # Check if file_paths is None or empty
    if file_paths is None or len(file_paths) == 0:
        print("No files found in the specified directory.")
    else:
        # Iterate over the files
        for file_path in file_paths:
            file_name = os.path.basename(file_path)
            parts = file_name.split('.')
            station_name = parts[1]
            longitude, latitude = get_station_coordinates(' ', full_coord_file_path, station_name)
        
            # Read MiniSEED file
            stream = read(file_path)
            print(file_name, stream)
        
            # Set trace start time and processing window
            reference_date = stream[0].stats.starttime
            processing_start_time = datetime(reference_date.year, reference_date.month, reference_date.day, 0, 0, 0)
            t1 = UTCDateTime(processing_start_time)
            t2 = t1 + 24 * 60 * 60
            stream.trim(starttime=t1, endtime=t2, pad=True, nearest_sample=False, fill_value=0)
        
            # Processing steps
            # 0. Decimate the signal to 10 Hz
            current_sampling_rate = stream[0].stats.sampling_rate
            decimation_factor = int(current_sampling_rate / 10)
            stream.decimate(decimation_factor)
        
            # Steps 1-6: Detrend, apply taper, band-pass filter, discard high-amplitude parts, and normalize
            for trace in stream:
                trace.detrend('demean')
                trace.detrend('linear')
                trace.taper(max_percentage=0.05, type='cosine')
                std_dev = np.std(trace.data)
                trace.data[np.abs(trace.data) > 3 * std_dev] = 0
                trace.data = np.sign(trace.data)
        
            # Step 7: Whiten the signal
            for trace in stream:
                trace = whiten_trace(trace, 0.1, 2.0)
                print('Whitening done!')
        
            # Step 8: Band-pass filtering again
            for trace in stream:
                trace.data = filtfilt(*butter(4, [0.1, 4.95], btype='band', fs=trace.stats.sampling_rate), trace.data)
        
            # Step 9: Remove instrument response using poles and zeros
            # for trace in stream:
            #     Make sure to replace the placeholders with actual instrument response information for your instruments
            #     paz_remove = {'poles': [1 + 2j, 1 - 2j], 'zeros': [0, 0], 'sensitivity': 1.0}
            #     trace.stats.response = {'response': AttribDict({'instrument_sensitivity': paz_remove})}
            #     trace.remove_response(output="VEL", pre_filt=[0.001, 0.005, 45, 50])
        
            # Write SAC file
            sac_trace = SACTrace.from_obspy_trace(stream[0])
            sac_trace.stla = float(latitude)
            sac_trace.stlo = float(longitude)
            component = "EHZ" if channels == "EHZ.D" else "EHZ.D"
            year = str(stream[0].stats.starttime.year)
            julian_day = '%03i' % (stream[0].stats.starttime.julday)
            hour = '%02i' % (stream[0].stats.starttime.hour)
            minute = '%02i' % (stream[0].stats.starttime.minute)
            second = '%02i' % (stream[0].stats.starttime.second)
            output_file_name = f'{station_name}.{year}.{julian_day}.{hour}.{minute}.{second}.{component}.sac'
            output_sac_path = os.path.join(input_sac_directory, output_file_name)
            sac_trace.write(output_sac_path)
            print(output_sac_path)

if __name__ == "__main__":
    input_mseed_directory = "/path/to/your/input_mseed_directory"
    input_sac_directory = "/path/to/your/input_sac_directory"
    full_coord_file_path = "/path/to/your/full_coord_file.txt"
    channels = ["EHZ.D"]

    process_mseed_files(input_mseed_directory, input_sac_directory, full_coord_file_path, channels)
