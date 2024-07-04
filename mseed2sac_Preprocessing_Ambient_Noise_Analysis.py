# miniseed2sac_Pre-Processing.py	
"""
Converts raw miniseed data into day-long files and converts them to SAC format.

Processing steps:
0. Remove Instrument Response (Optional)
1. Decimate the signal to 10 Hz
2. Remove the mean from each hour-long segment
3. Detrend the data
4. Apply a Tukey window to apply a 5% cosine Taper
5. Band-pass filtering between 0.10 and 5 Hz (0.2 - 10sec)
6. Discard signal parts with amplitude greater than 3 times their standard deviation
7. One-bit normalization
8. Whiten the signal
9. Band-pass filtering between 0.10 and 5 Hz (0.2 - 10sec)

Converts to SAC and adds station lat, lon info.

To get a left-handed system that matches IRIS, need:
- EHZ.K -> EHZ
"""

# Importing necessary libraries
import os
import numpy as np
from obspy import read, UTCDateTime
from obspy.io.sac import SACTrace
from datetime import datetime
from scipy.signal import butter, filtfilt
from obspy.signal.invsim import cosine_taper

# Function to whiten the input trace
def whiten(tr, freqmin, freqmax):
    """
    A cosine-tapered filter should be applied to the input trace in order to achieve whitening. 
    This can be achieved by spectral whitening of the trace 'tr' using a cosine tapered boxcar between 'freqmin' and 'freqmax'. 
    This method was developed by Gaia Soldati and Licia Faenza at the INGV (Soldati et al. (2015)

    Args:
        tr (obspy.Trace): The input trace.
        freqmin (float): The minimum frequency for the passband, in Hz.
        freqmax (float): The maximum frequency for the passband, in Hz.

    Returns:
        obspy.Trace: The whitened trace.
    """
    tr = tr.copy()
    npts = tr.stats.npts
    dt = tr.stats.delta
    nsmo = int(0.1 / dt)
    taper = cosine_taper(npts, p=0.05)
    tr.data *= taper
    FFTs = np.fft.fft(tr.data)
    freqs = np.fft.fftfreq(npts, d=dt)
    JJ = np.where((freqs >= freqmin) & (freqs <= freqmax))[0]
    FFTsW = np.zeros(npts, dtype=np.complex128)
    FFTsW[JJ[0]:JJ[-1]+1] = FFTs[JJ[0]:JJ[-1]+1]
    smo1 = (np.cos(np.linspace(np.pi/2., np.pi, nsmo+1))**2.).astype(np.float32)
    espo = np.exp(1j * np.angle(FFTs[JJ[0]:JJ[0]+nsmo+1]))
    FFTsW[JJ[0]:JJ[0]+nsmo+1] = smo1 * espo
    tr.data = np.fft.ifft(FFTsW).real
    tr.data /= np.abs(tr.data).max()
    return tr

# Function to get filenames based on file type
def get_filenames(sel, ext0, wd):
    """
    Get a list of filenames in a directory based on file type.

    Args:
        sel (int): Selection variable determining the file type.
        ext0 (str): The extension variable.
        wd (str): The selected path.

    Returns:
        list: List of filenames.
    """
    if sel == 1:  # 24bit CORE file
        ext1 = f'.r{ext0[1]}{ext0[2]}'
        ext2 = f'.R{ext0[1]}{ext0[2]}'
        allFilesPaths = []
        for root, dirs, files in os.walk(wd):
            for file in files:
                if file.endswith(ext1):
                    fls = os.path.join(root, file)
                    allFilesPaths.append(fls)
                if file.endswith(ext2):
                    fls = os.path.join(root, file)
                    allFilesPaths.append(fls)
    if sel != 1:  # Other file type
        ext1 = f'.{ext0}'
        allFilesPaths = []
        for root, dirs, files in os.walk(wd):
            for file in files:
                if file.endswith(ext1):
                    fls = os.path.join(root, file)
                    allFilesPaths.append(fls)
    return allFilesPaths

# Function to set station coordinates from a file
def sta_coord(delim, path, sta0):
    """
    Get station coordinates from a file.

    Args:
        delim (str): The delimiter used to split the line.
        path (str): The path of the STA_COORD file containing the necessary info.
        sta0 (str): The selected station ID.

    Returns:
        tuple: Latitude and Longitude of the station.
    """
    with open(path) as f:
        for line in f:
            l1 = line.strip().split('\t')  # Split the line using the tab character
            sta = l1[0]
            lat = l1[1]
            lon = l1[2]
            if sta == sta0:
                lat0 = lat
                lon0 = lon
                break
    return lat0, lon0

# Setup paths
path2mseed = "PATH/TO/MINISEED/FILES"  # Input mseed data path
path2sac = "PATH/TO/SAVE/SAC/FILES"  # Output sac path
path2sta = "PATH/TO/STATIONS/COORDINATES/INFORMATION"  # Station file (station ID, lat, lon, elevation)

chs = ["EHZ.K"]  # Channel names
pathlist = get_filenames(2, 'mseed', path2mseed)  # Get list of filenames

# Main function
if __name__ == "__main__":
    # Loop over all files
    for mpath in pathlist:
        mflnm = os.path.basename(mpath)
        var = mflnm.split('.')
        staN = var[1]
        lon, lat = sta_coord(' ', path2sta, staN)  # Get station coordinates
        st = read(mpath)  # Read mseed file
        tstart = st[0]
        refdate = st[0].stats.starttime
        tstart = datetime(refdate.year, refdate.month, refdate.day, 0, 0, 0)
        # Make sure trace is 24 hours long
        t1 = UTCDateTime(tstart)
        t2 = t1 + 24 * 60 * 60
        st.trim(starttime=t1, endtime=t2, pad=True, nearest_sample=False, fill_value=0)

        # Processing steps
        # 0. Remove Instrument Response
        #for tr in st:
            #paz_remove = {'gain': "instrument_gain", 'poles': "instrument_poles", 'zeros': "instrument_zeros", 'sensitivity': "instrument_sensitivity"}
            #trace.stats.response = {'response': {'instrument_sensitivity': paz_remove}}
            #trace.remove_response(output="VEL", pre_filt=[0.001, 0.005, 45, 50])
        current_sampling_rate = st[0].stats.sampling_rate
        decimation_factor = int(current_sampling_rate / 10)
        # 1. Decimate the signal to 10 Hz
        st.decimate(decimation_factor)
        # 2. Remove the mean from each hour-long segment
        st.detrend('demean')
        # 3. Detrend the data
        st.detrend('linear')
        # 4. Apply a Tukey window to apply a 5% cosine Taper
        st.taper(max_percentage=0.05, type='cosine')
        # 5. Band-pass filtering between 0.10 and 5 Hz (0.2 - 10sec)
        for tr in st:
            tr.data = filtfilt(*butter(4, [0.1, 4.95], btype='band', fs=tr.stats.sampling_rate), tr.data)
        # 6. Discard signal parts with amplitude greater than 3 times their standard deviation
        for tr in st:
            std = np.std(tr.data)
            tr.data[np.abs(tr.data) > 3 * std] = 0
        # 7. One-bit normalization
        for tr in st:
            tr.data = np.sign(tr.data)
        # 8. Whiten the signal
        for tr in st:
            tr = whiten(tr, 0.1, 2.0)
        # 9. Band-pass filtering between 0.10 and 5 Hz (0.2 - 10sec)
        for tr in st:
            tr.data = filtfilt(*butter(4, [0.1, 4.95], btype='band', fs=tr.stats.sampling_rate), tr.data)

        # Write SAC file
        sac = SACTrace.from_obspy_trace(st[0])
        sac.stla = float(lat)
        sac.stlo = float(lon)
        if chs == "EHZ.K":
            comp = "EHZ"
        else:
            comp = "EHZ.K"
        sac.kcmpnm = comp
        yr = str(st[0].stats.starttime.year)
        jday = '%03i' % (st[0].stats.starttime.julday)
        hr = '%02i' % (st[0].stats.starttime.hour)
        mn = '%02i' % (st[0].stats.starttime.minute)
        sec = '%02i' % (st[0].stats.starttime.second)
        pathout = f'{staN}.{yr}.{jday}.{hr}.{mn}.{sec}.{comp}.sac'
        sac_out = os.path.join(path2sac, pathout)
#-------#sac.id = f"{staN}.{comp}.."
        sac.write(sac_out)
        print(sac_out)
