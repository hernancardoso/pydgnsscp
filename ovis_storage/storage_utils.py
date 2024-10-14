# storage_utils.py
import binascii
from timestamps import timestamp_to_unix_millis
import json

def read_file(file_path):
    """
    Reads the RTCM message file and returns a list of dicts containing the UNIX timestamp, the message, and its decoded form.
    """
    with open(file_path, 'r') as file:
        reports = json.load(file)
        fixes = []
        for report in reports:
            try:
                timestamp = report['timestamp'] # Is alwasys present
                # la, lo, tg should be arrays of the same size
                tg = report.get("tg")
                la = report.get("la")
                lo = report.get("lo")
                alt = report.get("alt")
                sat = report.get("sat")

                if None in (tg, la, lo, alt, sat):
                    continue  # Skip this iteration

                sats = sat.split(",")

                for i in range(len(la)):
                    sats_data_used_for_fix = sats[i].split("/")
                    data = {}
                    for j in range(len(sats_data_used_for_fix)):
                        sat_data = sats_data_used_for_fix[j].split("-")
                        prn = sat_data[0]
                        az = int(sat_data[4])
                        el = int(sat_data[5])
                        data[int(prn)] = (az, el)

                    fix = {
                        "timestamp_unix": tg[i],
                        "timestamp": timestamp,
                        "la": la[i],
                        "lo": lo[i],
                        "alt": alt[i],
                        "sat": data
                    }
                    fixes.append(fix)
            except:
                continue

    return fixes

def save_message(file_path, message, timestamp):
    """
    Saves an RTCM message along with its timestamp to a file.
    """
    with open(file_path, 'a') as file:
        file.write(f"{timestamp};{message}\n")
