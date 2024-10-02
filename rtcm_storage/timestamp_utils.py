# timestamp_utils.py

import time
from datetime import datetime

def timestamp_to_unix(timestamp_str):
    """
    Converts a timestamp in the format 'YYYY-MM-DDTHH:MM:SS.ssssss' to UNIX time.
    """
    dt = datetime.strptime(timestamp_str, "%Y-%m-%dT%H:%M:%S.%f")
    return int(time.mktime(dt.timetuple()) + dt.microsecond / 1e6)


def find_closest_timestamp(target_timestamp, table):
    """
    Finds the closest timestamp in the table to the target timestamp.
    """
    return min(table, key=lambda x: abs(x['timestamp_unix'] - target_timestamp))
