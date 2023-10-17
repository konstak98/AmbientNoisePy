import os
from obspy import read, UTCDateTime
from obspy.io.sac import SACTrace
from obspy.core import AttribDict
from scipy.signal import butter, filtfilt
from obspy.signal.invsim import cosine_taper
import numpy as np
from datetime import datetime

# Function to whiten the input trace
def whiten(tr, freqmin, freqmax):
    # Implementation details omitted for brevity...
    pass

# Function to get filenames based on file type and extension
def get_filenames(sel, ext0, wd):
    if sel == 2:
        ext1 = f".{ext0}"
        allFilesPaths = [os.path.join(root, file) for root, dirs, files in os.walk(wd) for file in files if file.endswith(ext1)]
    else:
        ext1 = f".r{ext0[1]}{ext0[2]}"
        ext2 = f".R{ext0[1]}{ext0[2]}"
        allFilesPaths = [os.path.join(root, file) for root, dirs, files in os.walk(wd) for file in files if file.endswith(ext1) or file.endswith(ext2)]
    return allFilesPaths

# Function to retrieve station coordinates
def sta_coord(delim, path, sta0):
    with open(path) as f:
        for line in f:
            l1 = line.split(delim)
            sta = l1[0]
            lat = l1[1]
            lon = l1[2]
            if sta == sta0:
                return lat, lon

    # If station ID is not found, return None
    return None, None

# Setup paths
path2mseed = "/path/to/your/input_mseed_directory" # Replace these placeholders with your specific path
path2sac = "/path/to/your/input_sac_directory" # Replace these placeholders with your specific path
path2sta = "/path/to/your/full_coord_file.txt" # Replace these placeholders with your specific path
chs = ["EHZ.D"]
pathlist = get_filenames(2, 'mseed', path2mseed)

# Check if pathlist is None or empty
if pathlist is None or len(pathlist) == 0:
    print("No files found in the specified directory.")
else:
    # Iterate over the files
    for mpath in pathlist:
        mflnm = os.path.basename(mpath)
        var = mflnm.split('.')
        staN = var[1]
        lon, lat = sta_coord(' ', path2sta, staN)
    
        # Read mseed file
        st = read(mpath)
        print(mflnm, st)
    
        # Set trace start time and processing window
        refdate = st[0].stats.starttime
        tstart = datetime(refdate.year, refdate.month, refdate.day, 0, 0, 0)
        t1 = UTCDateTime(tstart)
        t2 = t1 + 24 * 60 * 60
        st.trim(starttime=t1, endtime=t2, pad=True, nearest_sample=False, fill_value=0)
    
        # Processing steps
        # 0. Decimate the signal to 10 Hz
        current_sampling_rate = st[0].stats.sampling_rate
        decimation_factor = int(current_sampling_rate / 10)
        st.decimate(decimation_factor)
    
        # Steps 1-6: Detrend, apply taper, band-pass filter, discard high-amplitude parts, and normalize
        for tr in st:
            tr.detrend('demean')
            tr.detrend('linear')
            tr.taper(max_percentage=0.05, type='cosine')
            std = np.std(tr.data)
            tr.data[np.abs(tr.data) > 3 * std] = 0
            tr.data = np.sign(tr.data)
    
        # Step 7: Whiten the signal
        for tr in st:
            tr = whiten(tr, 0.1, 2.0)
            print('Whitening done!')
    
        # Step 8: Band-pass filtering again
        for tr in st:
            tr.data = filtfilt(*butter(4, [0.1, 4.95], btype='band', fs=tr.stats.sampling_rate), tr.data)
    
        # Step 9: Remove instrument response using poles and zeros
        #for tr in st:
            # Make sure to replace the placeholders with actual instrument response information for your instruments
            #paz_remove = {'poles': [1 + 2j, 1 - 2j], 'zeros': [0, 0], 'sensitivity': 1.0}
            #tr.stats.response = {'response': AttribDict({'instrument_sensitivity': paz_remove})}
            #tr.remove_response(output="VEL", pre_filt=[0.001, 0.005, 45, 50])
    
        # Write SAC file
        sac = SACTrace.from_obspy_trace(st[0])
        sac.stla = float(lat)
        sac.stlo = float(lon)
        comp = "EHZ" if chs == "EHZ.D" else "EHZ.D"
        sac.kcmpnm = comp
        yr = str(st[0].stats.starttime.year)
        jday = '%03i' % (st[0].stats.starttime.julday)
        hr = '%02i' % (st[0].stats.starttime.hour)
        mn = '%02i' % (st[0].stats.starttime.minute)
        sec = '%02i' % (st[0].stats.starttime.second)
        pathout = f'{staN}.{yr}.{jday}.{hr}.{mn}.{sec}.{comp}.sac'
        sac_out = os.path.join(path2sac, pathout)
        sac.write(sac_out)
        print(sac_out)