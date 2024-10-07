# main.py
from dgnsscp.dgnsscp import DGNSSCP
from models import LatLonAltCoord
from rtcm_storage import read_file, count_rtcm_messages_by_id
from timestamps import timestamp_to_unix_millis, find_closest_timestamp


file_path = "timestamp.txt"

# Reading RTCM messages from the file (including decoding)
data_table = read_file(file_path)

# Output the data with the decoded message
for entry in data_table:
    print(f"Timestamp (UNIX): {entry['timestamp_unix']}")
    print(f"Raw RTCM Message: {entry['message']}")
    print("-" * 40)

test = count_rtcm_messages_by_id(data_table)
print(test)

#print("Starting PRC correction")
#dgnsscp = DGNSSCP(fixed_base_coords=LatLonAltCoord(lat=10, alt=10, lon=10))
#dgnsscp.process_messages(data_table)
#print(dgnsscp.prcs_history)
# Output the closest entry
#print(f"Closest timestamp entry: {closest_entry}")
