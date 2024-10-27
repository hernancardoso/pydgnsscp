# main.py
from dgnsscp.coordinate_transformations import ecef_to_geodetic
from dgnsscp.dgnsscp import DGNSSCP
from rtcm_storage import read_file, count_rtcm_messages_by_id
import ovis_storage
import folium


# Reading RTCM messages from the file (and also decode them)
data_table = read_file("timestamp.txt")

# Reading saved reports of collar
data = ovis_storage.read_file("reports.json")

# ECEF coordinates of UYMO mountpoint of REGNAROU
uymo_ecef = (2909132.9812000003, -4355451.2094, -3627801.3051) 
# Geodetic coordinates of UYMO
lat, lon, alt = ecef_to_geodetic(*uymo_ecef)

# Start DGNSSCP class
dgnsscp = DGNSSCP(fixed_base_coords = (lat, lon, alt))

# Iterate over RTCM messages generating the pseudo range correction.
# This step only uses data from the CORS, the idea is to use the UYMO as a rover.
# 1004 RTCM messages provide the observed pseudorange, and by knowing the exact position of the UYMO base station
# and the positions of the satellites it is possible to calculate the exact pseudorange, therefore 
# the substraction of these terms is the pseudo range correction for a particular satellite.
prcs = dgnsscp.process_messages(data_table)

# Now the rover is the collar, for each fix of the collar lets see which satellites it used, then
# apply the pseudo range correction to it and then recalculate the position with least square method
corrections = dgnsscp.correction_projection(data)

# lat,lon: are the original lat and lon calculated by the rover (the collar)
# lat_mod, lon_mod: are the fixes with the correction_projection algorithm 
for correction in corrections:
    print(f'Original lat: {correction['lat']}  - Fixed lat: {correction['lat_mod']}')
    print(f'Original lon: {correction['lon']}  - Fixed lon: {correction['lon_mod']}')
    print("-" * 30)