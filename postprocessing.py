# main.py
from dgnsscp.coordinate_transformations import ecef_to_geodetic
from dgnsscp.dgnsscp import DGNSSCP
from models import LatLonAltCoord
from rtcm_storage import read_file, count_rtcm_messages_by_id
import ovis_storage
import folium


# Reading RTCM messages from the file (including decoding)
data_table = read_file("timestamp.txt")


data = ovis_storage.read_file("reports.json")

(x, y, z) = (2909132.9812000003, -4355451.2094, -3627801.3051) # ECEF coordinates of UYMO
lat, lon, alt = ecef_to_geodetic(x, y, z)

fixed_base_geodetic = LatLonAltCoord(lat, lon, alt)

dgnsscp = DGNSSCP(fixed_base_coords = fixed_base_geodetic)
test = dgnsscp.process_messages(data_table)

corrections = dgnsscp.correction_projection(data)
print("End")

print(count_rtcm_messages_by_id(data_table))

latitudes = []
longitudes = []

latitudes_mod = []
longitudes_mod = []

for fix in corrections:
    la, lo = fix['lat'], fix['lon']
    la_mod, lo_mod = fix['lat_mod'], fix['lon_mod']
    latitudes.append(la)
    longitudes.append(lo)

    latitudes_mod.append(la_mod)
    longitudes_mod.append(lo_mod)


# Center the map on a location
m = folium.Map(location=[40.7128, -74.0060], zoom_start=2)  # Centered on New York for this example

for lat, lon in zip(latitudes, longitudes):
    folium.Marker(location=[lat, lon], icon=folium.Icon(color='red')).add_to(m)

for lat, lon in zip(latitudes, longitudes):
    folium.Marker(location=[lat, lon], icon=folium.Icon(color='red')).add_to(m)

for lat, lon in zip(latitudes_mod, longitudes_mod):
    folium.Marker(location=[lat, lon], icon=folium.Icon(color='green')).add_to(m)

# Display the map
m.save("map.html")

