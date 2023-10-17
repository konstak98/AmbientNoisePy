# AmbientNoisePy
Ambient-Noise Seismology Package for Analysis

```markdown
# mseed2sac_preprocessing_ambient_noise_analysis.py

This repository contains a specific Python script named `mseed2sac_preprocessing_ambient_noise_analysis.py` tailored for preprocessing seismic data in MiniSEED format and generating SAC files with specific naming conventions. Additional scripts for similar purposes may be added in the future.

## Table of Contents

- [Introduction](#introduction)
- [Requirements](#requirements)
- [Input MiniSEED Format Requirements](#input-miniseed-format-requirements)
- [Usage](#usage)
- [Functions](#functions)
- [License](#license)

## Introduction

The `mseed2sac_preprocessing_ambient_noise_analysis.py` script processes seismic data in MiniSEED format and applies a series of operations to each trace:
1. Decimate the signal to 10 Hz.
2. Detrend and apply a cosine taper.
3. Band-pass filter the data.
4. Discard high-amplitude parts.
5. Normalize the data.
6. Whiten the signal.
7. Band-pass filter the whitened signal.
8. Remove instrument response (commented out, requires actual instrument response information).
9. Write the processed data to a SAC file with a specific naming convention.

## Requirements

To use the script, the following dependencies are required:
- `obspy`: A Python framework for processing seismological data.
- `scipy`: A library for scientific and technical computing, including signal processing.
- `numpy`: A fundamental package for scientific computing with Python.

Install the dependencies using pip:

```bash
pip install obspy scipy numpy
```

## Input MiniSEED Format Requirements

The input MiniSEED files should adhere to the following naming convention:

```
SM.001.00.STATION_NAME.YEAR.JULIAN_DAY.mseed
```

- `SM`: Fixed prefix for the MiniSEED filename.
- `001`: Placeholder for a numbering scheme or any additional identifier.
- `00`: Placeholder for a numbering scheme or any additional identifier.
- `STATION_NAME`: Unique identifier for the station.
- `YEAR`: Year of the seismic data.
- `JULIAN_DAY`: Julian day of the year when the data was recorded.
- `mseed`: Fixed extension indicating the file format.

Ensure each MiniSEED file follows this naming convention for proper processing using the provided script.

## Usage

1. Adjust the paths and settings in the script to match your specific requirements, including the desired MiniSEED filename format.
2. Ensure the required dependencies are installed (see Requirements).
3. Execute the script using Python:

```bash
python mseed2sac_preprocessing_ambient_noise_analysis.py
```

## Functions

### `whiten(tr, freqmin, freqmax)`

This function whitens the input trace using specified frequency bounds.

### `get_filenames(sel, ext0, wd)`

This function retrieves filenames based on the file type and extension.

### `sta_coord(delim, path, sta0)`

This function retrieves station coordinates based on a specified station ID.

## License

This script is provided under the [MIT License](LICENSE). Feel free to modify and distribute it as needed.
```
