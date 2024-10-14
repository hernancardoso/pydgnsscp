# timestamp_utils.py

import time
from datetime import datetime

def timestamp_to_unix_millis(timestamp_str):
    """
    Converts a timestamp in the format 'YYYY-MM-DDTHH:MM:SS.ssssss' to UNIX time, keeping only milliseconds.
    """
    dt = datetime.strptime(timestamp_str, "%Y-%m-%dT%H:%M:%S.%f")
    unix_timestamp = int(dt.timestamp() * 1000)  # Convert to milliseconds and cast to int
    return unix_timestamp

def find_closest_timestamp(target_timestamp, table):
    """
    Finds the closest timestamp in the table to the target timestamp.
    The search is going to be indexed by the 'timestamp_unix' column
    """
    return min(table, key=lambda x: abs(x['timestamp_unix'] - target_timestamp))

def unix_timestamp_millis_to_datetime(unix_timestamp_millis):
    return datetime.fromtimestamp(unix_timestamp_millis / 1000)