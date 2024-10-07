import math

from dgnsscp.coordinate_transformations import calculate_enu_from_geodetic_ecef, ecef_to_geodetic, geodetic_to_ecef, \
    calculate_ecef_satellite_vector
from dgnsscp.dgnsscp_utils import retrieve_ephemeris, calculate_ecef_satellite_position, calculate_azimuth_elevation
from models import LatLonAltCoord, AzElCoord


class DGNSSCP:
    def __init__(self, fixed_base_coords: LatLonAltCoord):
        self.fixed_base_coords = fixed_base_coords
        self.prcs_history = []  # Array of objects { timestamp_unix, prcs }


    def process_messages(self, messages):
        np_fixed_base_coord = self.fixed_base_coords.to_numpy()
        satellites_az_ez = {}

        for message in messages:
            try:
                timestamp = message["timestamp_unix"]
                rtcm_message = message["message"]

                if rtcm_message:
                    if rtcm_message.identity == '1019':
                        satellite_id = rtcm_message.DF009
                        ephemeris = retrieve_ephemeris(rtcm_message)

                        np_satellite_ecef = calculate_ecef_satellite_position(ephemeris)

                        enu = calculate_enu_from_geodetic_ecef(np_fixed_base_coord, np_satellite_ecef)
                        azimuth, elevation = calculate_azimuth_elevation(enu)

                        satellites_az_ez[satellite_id] = AzElCoord(azimuth, elevation)

                        print(f"azimuth: {azimuth}, elevation: {elevation}")

                    if rtcm_message.identity == '1005' or rtcm_message.identity == '1006':
                        x = rtcm_message.DF025 * 0.0001
                        y = rtcm_message.DF026 * 0.0001
                        z = rtcm_message.DF027 * 0.0001

                        lat, lon, alt = ecef_to_geodetic(x, y, z)

                        observed_base_coords = LatLonAltCoord(lat, lon, alt)

                        prcs = {}

                        for prn, az_el in satellites_az_ez.items():  # used for the fix
                            prc = generate_prc_for_satellite(observed_base_coords, self.fixed_base_coords, az_el.azimuth, az_el.elevation, 0)
                            print("El PRC es :", prc)
                            prcs[prn] = prc

                        row = {
                            "timestamp_unix": timestamp,
                            "prcs": prcs
                        }
                        self.prcs_history.append(row)
                        print(f"Latitude: {lat}, Longitude: {lon}, Altitude: {alt}")
                    # else:
                    #    print("Message does not contain ECEF coordinates.")
                else:
                    print("Failed to parse RTK message.")
            except Exception as e:
                print(f"Error parsing RTK message: {e}")


def generate_prc_for_satellite(observed_gga: LatLonAltCoord, base_coord: LatLonAltCoord, sat_acimut, sat_elevacion,
                                sat_residuals=0):
    # def generate_prc_for_satellite(observed_gga: GPGGA, base_station_pos, gsv):

    # Obtengo las coordenadas ECEF en funcion del mensaje NMEA GPGGA.
    x_medido, y_medido, z_medido = geodetic_to_ecef(observed_gga.lat, observed_gga.lon,
                                                                    observed_gga.alt)

    # Obtengo las coordenadas ECEF en funcion de la latitud, longitud y altitud conocidas de la base station.
    x_conocido, y_conocido, z_conocido = geodetic_to_ecef(base_coord.lat, base_coord.lon,
                                                                          base_coord.alt)

    # Obtengo las coordenadas ECEF del satelite i "ideal", usando las coordenadas de la base station
    x_sat_ideal, y_sat_ideal, z_sat_ideal = calculate_ecef_satellite_vector(base_coord.lat,
                                                                                            base_coord.lon,
                                                                                            base_coord.alt, sat_acimut,
                                                                                            sat_elevacion)

    # Obtengo las coordenadas ECEF del satelite i, usando las coordenadas calculadas para la base station
    x_sat, y_sat, z_sat = calculate_ecef_satellite_vector(observed_gga.lat, observed_gga.lon,
                                                                          observed_gga.alt, sat_acimut, sat_elevacion)

    # Con las ECEF del satelite y las ECEF del receptor puedo calcular el pseudorango ideal y el observado.

    p_ideal = math.sqrt(
        (x_sat_ideal - x_conocido) ** 2 + (y_sat_ideal - y_conocido) ** 2 + (z_sat_ideal - z_conocido) ** 2)

    p_observado = math.sqrt((x_sat - x_medido) ** 2 + (y_sat - y_medido) ** 2 + (z_sat - z_medido) ** 2)

    # Para aplicar los residuals los sumo al pseudorango observado.
    # Primero tengo que identificar el correcto residual.
    # El orden de los PRN de arr_sat_gsv es el mismo que el de arr_residuals.

    p_observado = p_observado + sat_residuals
    prc = p_ideal - p_observado

    return prc
