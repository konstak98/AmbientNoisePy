import os
import numpy as np
from obspy import read
from obspy.signal.cross_correlation import correlate

# Define constants
# Insert your data directory path here
data_dir = "INSERT_YOUR_DATA_DIRECTORY_PATH"
window = 1
maxlag = 1500
lags = np.arange(-maxlag, maxlag, window)

def compute_cross_correlations(data_dir, day0, dayn):
    # Get a list of all SAC files in the data directory
    sac_files = sorted([f for f in os.listdir(data_dir) if f.endswith(".sac") and "EHZ.D" in f])
    print(sac_files)

    it = 1
    ndays = dayn - day0 + 1
    nsta = 80

    # Get unique station IDs
    sta_id_list = list(set([fl[:3] for fl in sac_files]))

    for sta1 in sta_id_list:
        sta1_int = int(sta1)

        for day in range(ndays):
            d0 = day0 + day
            true_day = d0
            flnm0 = '{}.2013.{:03d}.00.00.00.EHZ.D.sac'.format(sta1, d0)

            for sta2 in sta_id_list:
                sta2_int = int(sta2)
                if sta2_int < sta1_int or sta1 == sta2:
                    continue

                flnm1 = '{}.2013.{:03d}.00.00.00.EHZ.D.sac'.format(sta2, d0)

                filename1 = os.path.join(data_dir, flnm0)
                filename2 = os.path.join(data_dir, flnm1)

                if not os.path.exists(filename2):
                    print('*** STATION {} for DAY {:03d} does not exist.'.format(sta2, d0))
                    continue

                st1 = read(filename1)
                st2 = read(filename2)
                tr1 = st1[0]
                tr2 = st2[0]

                # Calculate cross-correlation
                rs = correlate(tr1.data, tr2.data, maxlag, demean=True, normalize='naive', method='auto')

                # Adjust lags and rs to have the same length
                lags = np.arange(-maxlag, maxlag)
                rs = rs[:len(lags)]

                # Save cross-correlation to a file
                # Insert your output directory path here
                filename_ccf = "INSERT_YOUR_OUTPUT_DIRECTORY_PATH/crosscorrelation_{sta1}_{sta2}_{true_day:3d}.txt"
                np.savetxt(filename_ccf, rs)

def main():
    day0 = 248
    dayn = 248
    compute_cross_correlations(data_dir, day0, dayn)

if __name__ == "__main__":
    main()
