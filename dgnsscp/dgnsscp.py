import math
import numpy as np
from datetime import datetime

from gps.az_el_utils import SatellitePositionCalculator
from timestamps import unix_timestamp_millis_to_datetime, find_closest_timestamp
from dgnsscp.coordinate_transformations import geodetic_to_ecef, ecef_to_geodetic
from models import LatLonAltCoord


class DGNSSCP:
    def __init__(self, fixed_base_coords: LatLonAltCoord):
        """
        Initialize the DGNSSCP class with fixed base coordinates and navigation file.

        :param fixed_base_coords: LatLonAltCoord object with base station's latitude, longitude, and altitude.
        :param nav_file_path: Path to the navigation file used by the SatellitePositionCalculator.
        """
        self.fixed_base_coords = fixed_base_coords
        self.prcs_history = []  # List of dictionaries with keys 'timestamp_unix' and 'prcs'
        self.satpos_calculator = SatellitePositionCalculator("./gps/brdc2680.24n")


    def process_messages(self, messages):
        """
        Process a list of RTCM messages to calculate pseudorange corrections (PRCs).

        :param messages: List of dictionaries containing 'timestamp_unix' and 'message' (RTCM message).
        :return: List of PRC history entries.
        """
        for message in messages:
            try:
                timestamp = message["timestamp_unix"]
                datetime_timestamp = unix_timestamp_millis_to_datetime(timestamp)
                rtcm_message = message["message"]

                if rtcm_message:
                    if rtcm_message.identity == '1004':
                        prcs = self.calculate_prcs(rtcm_message, datetime_timestamp)
                        self.prcs_history.append({
                            'timestamp_unix': timestamp,
                            'prcs': prcs
                        })
                else:
                    print("Failed to parse RTCM message.")
            except Exception as e:
                print(f"Error parsing RTCM message: {e}")

        return self.prcs_history

    def correction_projection(self, reports):
        """
        Apply corrections to a list of reports to compute modified coordinates.

        :param reports: List of dictionaries containing report data from devices.
        """
        corrected_gnss = []
        for report in reports:
            satellite_data = []
            num_satellites = len(report["sat"])
            if(num_satellites < 4):
                continue # Not enough satellites for doing the fix
            timestamp = report["timestamp_unix"]

            # Measured coordinates from the collar
            lat = report["la"]
            lon = report["lo"]
            alt = report["alt"]

            # Convert geodetic to ECEF coordinates
            x_measured, y_measured, z_measured = geodetic_to_ecef(lat, lon, alt)

            # Get satellite ECEF positions at the given timestamp
            sat_ecef = self.satpos_calculator.calculate_all_positions(
                unix_timestamp_millis_to_datetime(timestamp)
            )

            # Find the closest PRC data to the report timestamp
            prc_data = find_closest_timestamp(timestamp, self.prcs_history)

            # For each satellite used in the fix, reconstruct the observed pseudorange and apply correction
            for prn in report["sat"]:
                try:
                    x_sat, y_sat, z_sat = sat_ecef[prn]
                except KeyError:
                    print(f"Satellite PRN {prn} position not found.")
                    continue

                # Calculate observed pseudorange
                p_observed = math.sqrt(
                    (x_sat - x_measured) ** 2 +
                    (y_sat - y_measured) ** 2 +
                    (z_sat - z_measured) ** 2
                )

                # Apply correction if available
                try:
                    prc = prc_data["prcs"][prn]
                    p_corrected = p_observed + prc
                except KeyError:
                    print(f"PRC for PRN {prn} not found.")
                    p_corrected = p_observed

                satellite_data.append((x_sat, y_sat, z_sat, p_corrected))

            # Now, given the satellites positions and fixed pseudorange 
            # do least squares to calculate the lateration
        
            x_mod, y_mod, z_mod = x_measured, y_measured, z_measured
            tau_mod = 0  # Initial clock bias
            H = np.zeros((num_satellites, 4))
            DeltaP = np.zeros(num_satellites)

            # Iterative least squares solution
            for _ in range(20):
                for j, (x_sat, y_sat, z_sat, pseudorange) in enumerate(satellite_data):
                    p = math.sqrt(
                        (x_sat - x_mod) ** 2 +
                        (y_sat - y_mod) ** 2 +
                        (z_sat - z_mod) ** 2
                    ) + tau_mod

                    #H, which is the design matrix
                    # stores partial derivatives of the pseudorange equations with respect 
                    # to the unknowns (x, y, z)
                    H[j] = np.array([
                        (x_mod - x_sat) / p,
                        (y_mod - y_sat) / p,
                        (z_mod - z_sat) / p,
                        299792458  # Speed of light in m/s
                    ]) 

                    DeltaP[j] = pseudorange - p #  store the differences between measured and predicted pseudoranges

                # Compute correction using least squares
                HtH = H.T @ H
                HtH_inv = np.linalg.pinv(HtH) # Pseudo inverse is being used bc the HtH matrix is poorly conditioned
                #HtH_inv = np.linalg.pinv(HtH) # Pseudo inverse is being used bc the HtH matrix is poorly conditioned

                HtDeltaP = H.T @ DeltaP
                DeltaX = HtH_inv @ HtDeltaP

                # Update estimates
                x_mod += DeltaX[0]
                y_mod += DeltaX[1]
                z_mod += DeltaX[2]
                tau_mod += DeltaX[3]

            # Convert ECEF to geodetic coordinates
            lat_mod, lon_mod, alt_mod = ecef_to_geodetic(x_mod, y_mod, z_mod)

            # Update report with modified coordinates
            report.update({
                "x_mod": x_mod,
                "y_mod": y_mod,
                "z_mod": z_mod,
                "lat_mod": lat_mod,
                "lon_mod": lon_mod,
                "alt_mod": alt_mod
            })
            corrected_gnss.append({
                "lat": lat,
                "lon": lon,
                "lat_mod": lat_mod,
                "lon_mod": lon_mod
            })
        return 


    def calculate_prcs(self, rtcm_message, datetime_timestamp):
        """
        Calculate pseudorange corrections from an RTCM message.

        :param rtcm_message: RTCM message containing satellite observations.
        :param datetime_timestamp: Timestamp of the message as a datetime object.
        :return: Dictionary of PRCs with PRN as keys.
        """
        sat_ecef = self.satpos_calculator.calculate_all_positions(datetime_timestamp)
        prcs = {}

        for i in range(1, 32):  # PRNs from 1 to 31
            pseudorange_field = f'DF011_{i:02d}'
            prn_field = f'DF009_{i:02d}'

            try:
                prn = getattr(rtcm_message, prn_field)
                pseudorange_raw = getattr(rtcm_message, pseudorange_field)

                # Estimate total range
                approximate_total_range = self.estimate_total_range(
                    sat_ecef, prn
                )

                # Compute number of whole milliseconds
                N = int(round((approximate_total_range - pseudorange_raw) / 299792.458))

                # Compute full pseudorange
                pseudorange_full = pseudorange_raw + N * 299792.458  # In meters

                # Generate PRC
                correction = self.generate_prc_for_satellite(
                    pseudorange_full, sat_ecef, prn
                )
                prcs[prn] = correction

            except AttributeError:
                # Attribute does not exist, move to next PRN
                continue
            except KeyError:
                print(f"PRN {prn} data not found in satellite positions.")
                continue
            except Exception as e:
                print(f"Error processing PRN {prn}: {e}")
                continue

        return prcs

    def estimate_total_range(self, sat_ecef, prn):
        """
        Estimate the total range from the base station to the satellite.

        :param sat_ecef: Dictionary of satellite ECEF coordinates.
        :param prn: PRN number of the satellite.
        :return: Estimated total range in meters.
        """
        # Convert base station geodetic coordinates to ECEF
        x_base, y_base, z_base = geodetic_to_ecef(
            self.fixed_base_coords.lat,
            self.fixed_base_coords.lon,
            self.fixed_base_coords.alt
        )

        base_ecef = np.array([x_base, y_base, z_base])
        satellite_ecef = np.array(sat_ecef[prn])

        # Compute Euclidean distance
        approximate_total_range = np.linalg.norm(satellite_ecef - base_ecef)

        return approximate_total_range

    def generate_prc_for_satellite(self, pseudorange_full, sat_ecef, prn):
        """
        Generate pseudorange correction for a satellite.

        :param pseudorange_full: Full pseudorange measurement.
        :param sat_ecef: Dictionary of satellite ECEF coordinates.
        :param prn: PRN number of the satellite.
        :return: Pseudorange correction value.
        """
        # Known base station ECEF coordinates
        x_known, y_known, z_known = geodetic_to_ecef(
            self.fixed_base_coords.lat,
            self.fixed_base_coords.lon,
            self.fixed_base_coords.alt
        )

        x_sat, y_sat, z_sat = sat_ecef[prn]

        # Ideal pseudorange based on known positions
        p_ideal = math.sqrt(
            (x_sat - x_known) ** 2 +
            (y_sat - y_known) ** 2 +
            (z_sat - z_known) ** 2
        )

        # Observed pseudorange
        p_observed = pseudorange_full

        # Compute pseudorange correction
        prc = p_ideal - p_observed

        return prc
