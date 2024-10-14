from pymap3d import ecef2geodetic, geodetic2ecef, ecef2enu
import math

def ecef_to_geodetic(x, y, z):
    lat, lon, alt = ecef2geodetic(x, y, z)
    return lat, lon, alt

def geodetic_to_ecef(lat, lon, alt):
    x, y, z = geodetic2ecef(lat = lat, lon = lon, alt = alt)
    return x, y, z

def ecef_to_enu(x, y, z, lat, lon, alt):
    East, North, Up = ecef2enu(x, y, z, lat, lon, alt)
    return East, North, Up

def ecef_satelite(acimut, elevacion, lat, lon, alt):  # listo
    """
    Calcula las coordenadas ECEF de un satélite en función del acimut, elevación, y la posición del receptor.

    Parámetros:
    lat -- Latitud del receptor en grados.
    lon -- Longitud del receptor en grados.
    alt -- Altitud del receptor sobre el nivel del mar en metros.
    azimuth -- Acimut del satélite en grados.
    elevation -- Elevación del satélite en grados.

    Asumo que la distancia del receptor al satéltie es 20200000 m.

    Devuelve:
    Una tupla con las coordenadas (x, y, z) del satélite en metros en el sistema ECEF.
    """

    distance = 20200000

    # Convertir acimut y elevación a radianes
    acimut_rad = math.radians(acimut)
    elevacion_rad = math.radians(elevacion)

    # Convertir latitud y longitud del receptor a radianes
    lat_rad = math.radians(lat)
    lon_rad = math.radians(lon)

    # Convertir las coordenadas del receptor a ECEF
    x_r, y_r, z_r = geodetic_to_ecef(lat, lon, alt)

    # Calcular los componentes del vector del satélite en el sistema de coordenadas del receptor
    x_s = distance * math.cos(elevacion_rad) * math.sin(acimut_rad)
    y_s = distance * math.cos(elevacion_rad) * math.cos(acimut_rad)
    z_s = distance * math.sin(elevacion_rad)

    # Matriz de rotación del sistema de coordenadas del receptor a ECEF
    R = [
        [-math.sin(lon_rad), math.cos(lon_rad), 0],
        [-math.sin(lat_rad) * math.cos(lon_rad), -math.sin(lat_rad) * math.sin(lon_rad), math.cos(lat_rad)],
        [math.cos(lat_rad) * math.cos(lon_rad), math.cos(lat_rad) * math.sin(lon_rad), math.sin(lat_rad)]
    ]

    # Transformar las coordenadas del satélite al sistema ECEF
    x_s_ecef = x_r + R[0][0] * x_s + R[0][1] * y_s + R[0][2] * z_s
    y_s_ecef = y_r + R[1][0] * x_s + R[1][1] * y_s + R[1][2] * z_s
    z_s_ecef = z_r + R[2][0] * x_s + R[2][1] * y_s + R[2][2] * z_s

    return x_s_ecef, y_s_ecef, z_s_ecef
