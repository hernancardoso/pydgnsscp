# main.py

from rtcm_storage import read_file, timestamp_to_unix, find_closest_timestamp

file_path = "data.txt"

# Reading RTCM messages from the file (including decoding)
data_table = read_file(file_path)

# Output the data with the decoded message
for entry in data_table:
    print(f"Timestamp (UNIX): {entry['timestamp_unix']}")
    print(f"Raw RTCM Message: {entry['message']}")
    print("-" * 40)

# Example: Find the closest timestamp to a given one
target_timestamp = timestamp_to_unix("2024-09-24T17:50:00.000000")
closest_entry = find_closest_timestamp(target_timestamp, data_table)

# Output the closest entry
print(f"Closest timestamp entry: {closest_entry}")
