import math
import random
from datetime import datetime
from timestamps import unix_timestamp_millis_to_datetime, find_closest_timestamp
import gps
import numpy as np
from dgnsscp.coordinate_transformations import geodetic_to_ecef, ecef_satelite
from models import LatLonAltCoord, AzElCoord


class DGNSSCP:
    def __init__(self, fixed_base_coords: LatLonAltCoord):
        self.fixed_base_coords = fixed_base_coords
        self.prcs_history = []  # Array of objects { timestamp_unix, prcs }
        self.observed_base_pseudoranges_history_list = []

    def process_messages(self, messages):
        visited = False
        az_el_dict = {}

        for message in messages:
            try:
                timestamp = message["timestamp_unix"]
                rtcm_message = message["message"]

                if rtcm_message:
                    # Message type for base station position,
                    # From this messages the az/el can be calculated
                    if rtcm_message.identity == '1005' or rtcm_message.identity == '1006':
                        if visited:
                            continue
                        else:
                            visited = True
                        x = rtcm_message.DF025
                        y = rtcm_message.DF026
                        z = rtcm_message.DF027
                        print(x,y,z)
                        datetime_timestamp = unix_timestamp_millis_to_datetime(timestamp)
                        navfile = "/Users/hernancardoso/pydgnsscp/gps/brdc2680.24n"
                        az_el_dict = gps.calculate_az_el_dict(receiver_pos_llh=(x, y, z), obs_time= datetime_timestamp, navfile = navfile)

                    if rtcm_message.identity == '1004': #PRC Received from satellite
                        # save the PRN and pseudorange in a dict
                        prcs = {}

                        # Assuming there are up to 32 satellites in this message (based on the field names you provided)
                        for i in range(1, 32):  # Loop through DF009_01 to DF009_09 for PRNs
                            # Construct the field names dynamically
                            prn_field = f'DF009_0{i}'
                            pseudorange_field = f'DF011_0{i}'

                            # Extract the satellite PRN and pseudorange (assuming DF011 is the pseudorange field)
                            try:
                                prn = getattr(rtcm_message, prn_field)
                                pseudorange = getattr(rtcm_message, pseudorange_field)

                                (az, el) = az_el_dict[prn]

                                correction = generate_prc_for_satellite(pseudorange, self.fixed_base_coords, az, el)
                                # Save the PRN and pseudorange in the dictionary
                                prcs[prn] = correction
                            except:
                                break

                        self.observed_base_pseudoranges_history_list.append({
                            'timestamp_unix': timestamp,
                            'prcs': prcs
                        })
                else:
                    print("Failed to parse RTK message.")
            except Exception as e:
                print(f"Error parsing RTK message: {e}")

        return self.observed_base_pseudoranges_history_list

    def correction_projection(self, reports):
        parsed_reports = reports
        x_fixed, y_fixed, z_r = geodetic_to_ecef(lat = self.fixed_base_coords.lat, lon = self.fixed_base_coords.lon, alt = self.fixed_base_coords.alt)

        for report in reports:
            arr_nodo = []

            number_of_satellites_used = len(report["sat"])
            timestamp = report["timestamp_unix"]

            lat = report["la"]
            lon = report["lo"]
            alt = report["alt"]
            x_medido, y_medido, z_medido = geodetic_to_ecef(lat, lon, alt)
            data = find_closest_timestamp(timestamp, self.observed_base_pseudoranges_history_list)
            for prn in report["sat"]:
                az, el = report["sat"][prn]
                x_sat, y_sat, z_sat = ecef_satelite(az, el, lat,lon,alt)
                p_observado = math.sqrt(( x_sat-x_medido)**2 + (y_sat-y_medido)**2 + (z_sat-z_medido)**2 )
                p_corregido = p_observado +  data["prcs"][prn] - random.randint(10000, 1000000)
                arr_nodo.append(( x_sat, y_sat, z_sat, p_corregido))

            x_mod, y_mod, z_mod = geodetic_to_ecef(lat, lon, alt)
            # Suponemos que el DeltaTau = 0 en el comienzo.
            tau_mod = 0
            H = np.zeros((number_of_satellites_used, 4))
            DeltaP = np.zeros(number_of_satellites_used)

            for h in range(20):
                # En este for armamos la matriz H fila por fila.
                for j in range(number_of_satellites_used):
                    x, y, z, pseudorango = arr_nodo[j]
                    p = math.sqrt(
                        (x - x_mod) ** 2 + (y - y_mod) ** 2 + (z - z_mod) ** 2) + tau_mod
                    H[j] = np.array(
                        [(x_mod - x) / p, (y_mod - y) / p, (z_mod - z) / p, 299792458]) # 299792458 speed light

                    # DeltaP = p_corregido - p_modelado
                    DeltaP[j] = pseudorango - p

                    # Calcular H^tH
                HtH = np.dot(H.T, H)

                # Calcular la inversa de H^tH
                HtH_inv = np.linalg.inv(HtH)

                # Calcular H^tP
                HtDeltaP = np.dot(H.T, DeltaP)

                # Calcular X
                DeltaX = np.dot(HtH_inv, HtDeltaP)

                x_mod = x_mod + DeltaX[0]
                y_mod = y_mod + DeltaX[1]
                z_mod = z_mod + DeltaX[2]
                tau_mod = tau_mod + DeltaX[3]

            #return x_mod, y_mod, z_mod, tau_mod
            report["x_mod"] = x_mod
            report["y_mod"] = y_mod
            report["z_mod"] = z_mod

def generate_prc_for_satellite(p_observado, base_coord: LatLonAltCoord, sat_acimut, sat_elevacion,
                                sat_residuals=0):
    # Obtengo las coordenadas ECEF en funcion de la latitud, longitud y altitud conocidas de la base station.
    x_conocido, y_conocido, z_conocido = geodetic_to_ecef(base_coord.lat, base_coord.lon,
                                                                          base_coord.alt)

    # Obtengo las coordenadas ECEF del satelite i "ideal", usando las coordenadas de la base station
    x_sat_ideal, y_sat_ideal, z_sat_ideal = ecef_satelite(sat_acimut, sat_elevacion, base_coord.lat, base_coord.lon, base_coord.alt)


    # Con las ECEF del satelite y las ECEF del receptor puedo calcular el pseudorango ideal y el observado.

    p_ideal = math.sqrt(
        (x_sat_ideal - x_conocido) ** 2 + (y_sat_ideal - y_conocido) ** 2 + (z_sat_ideal - z_conocido) ** 2)


    # Para aplicar los residuals los sumo al pseudorango observado.
    # Primero tengo que identificar el correcto residual.
    # El orden de los PRN de arr_sat_gsv es el mismo que el de arr_residuals.

    p_observado = p_observado + sat_residuals
    prc = p_ideal - p_observado

    return prc


