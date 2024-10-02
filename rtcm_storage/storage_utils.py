# storage_utils.py
from rtcm_storage import timestamp_to_unix
from rtcm_storage.processing_utils import decode_rtcm_message


def read_file(file_path):
    """
    Reads the RTCM message file and returns a list of dicts containing the UNIX timestamp, the message, and its decoded form.
    """
    table = []
    with open(file_path, 'r') as file:
        for line in file:
            timestamp_str, message = line.strip().split(';')
            unix_timestamp = timestamp_to_unix(timestamp_str)
            decoded_message = decode_rtcm_message(message)
            table.append({
                "timestamp_unix": unix_timestamp,
                "message": message,
                "decoded_message": decoded_message
            })
    return table

def save_message(file_path, message, timestamp):
    """
    Saves an RTCM message along with its timestamp to a file.
    """
    with open(file_path, 'a') as file:
        file.write(f"{timestamp};{message}\n")
