# main.py
from dgnsscp.coordinate_transformations import ecef_to_geodetic, geodetic_to_ecef
from dgnsscp.dgnsscp import DGNSSCP
from models import LatLonAltCoord
from rtcm_storage import read_file, count_rtcm_messages_by_id
import ovis_storage
from timestamps import timestamp_to_unix_millis, find_closest_timestamp

file_path = "timestamp.txt"

# Reading RTCM messages from the file (including decoding)
data_table = read_file(file_path)


data = ovis_storage.read_file("/Users/hernancardoso/pydgnsscp/reports.json")

(x, y, z) = (2909132.9812000003, -4355451.2094, -3627801.3051)
lat, lon, alt = ecef_to_geodetic(x, y, z)

fixed_base_geodetic = LatLonAltCoord(lat, lon, alt)


dgnsscp = DGNSSCP(fixed_base_coords = fixed_base_geodetic)
test = dgnsscp.process_messages(data_table)

corrections = dgnsscp.correction_projection(data)
print("End")

print(count_rtcm_messages_by_id(data_table))
#print(dgnsscp.prcs_history)
# Output the closest entry
#print(f"Closest timestamp entry: {closest_entry}")
