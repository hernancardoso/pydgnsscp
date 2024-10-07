# storage_utils.py
import binascii
from .processing_utils import decode_rtcm_message
from timestamps import timestamp_to_unix_millis


def read_file(file_path):
    """
    Reads the RTCM message file and returns a list of dicts containing the UNIX timestamp, the message, and its decoded form.
    """
    table = []
    with open(file_path, 'r') as file:
        for line in file:
            timestamp_str, message_hex = line.strip().split(';')
            byte_message = binascii.unhexlify(message_hex)
            unix_timestamp = timestamp_to_unix_millis(timestamp_str)

            message_id, decoded_message = decode_rtcm_message(byte_message)
            table.append({
                "timestamp_unix": unix_timestamp,
                "byte_message": byte_message,
                "message_id": message_id,
                "message": decoded_message
            })
    return table

def save_message(file_path, message, timestamp):
    """
    Saves an RTCM message along with its timestamp to a file.
    """
    with open(file_path, 'a') as file:
        file.write(f"{timestamp};{message}\n")
