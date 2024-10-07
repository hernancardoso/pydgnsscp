import math
from typing import Tuple
import numpy as np

a = 6378137.0  # Semieje mayor en metros
f = 1 / 298.257223563  # Aplanamiento
e2 = 2 * f - f ** 2  # Excentricidad al cuadrado
eps = 1e-12  # Tolerancia para la convergencia de la iteración de Newton
distance = 20200000


def r_matrix(lon_rad, lat_rad):
    R = np.array([
        [-math.sin(lon_rad), math.cos(lon_rad), 0],
        [-math.sin(lat_rad) * math.cos(lon_rad), -math.sin(lat_rad) * math.sin(lon_rad), math.cos(lat_rad)],
        [math.cos(lat_rad) * math.cos(lon_rad), math.cos(lat_rad) * math.sin(lon_rad), math.sin(lat_rad)]
    ])
    return R


def calculate_ecef_satellite_vector(lat, lon, alt, acimut, elevacion):
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

    R = r_matrix(lon_rad, lat_rad)

    # Transformar las coordenadas del satélite al sistema ECEF
    x_s_ecef = x_r + R[0][0] * x_s + R[0][1] * y_s + R[0][2] * z_s
    y_s_ecef = y_r + R[1][0] * x_s + R[1][1] * y_s + R[1][2] * z_s
    z_s_ecef = z_r + R[2][0] * x_s + R[2][1] * y_s + R[2][2] * z_s

    return x_s_ecef, y_s_ecef, z_s_ecef


def geodetic_to_ecef(lat, lon, alt) -> np.ndarray:
    """
    Convierte coordenadas geodésicas (latitud, longitud, altitud MSL) a coordenadas ECEF (x, y, z).

    Devuelve:
    Una tupla con las coordenadas (x, y, z) en metros.
    """
    # Convertir latitud y longitud de grados a radianes
    lat_rad = math.radians(lat)
    lon_rad = math.radians(lon)

    # Radio de curvatura en el primer vertical
    N = a / math.sqrt(1 - e2 * math.sin(lat_rad) ** 2)

    # Coordenadas ECEF
    x = (N + alt) * math.cos(lat_rad) * math.cos(lon_rad)
    y = (N + alt) * math.cos(lat_rad) * math.sin(lon_rad)
    z = (N * (1 - e2) + alt) * math.sin(lat_rad)

    return np.array([x, y, z])


def ecef_to_geodetic(x, y, z):
    """
    Convierte coordenadas ECEF (x, y, z) a coordenadas geodésicas (latitud, longitud, altitud).
    x, y, z en metros.

    Tiene un pequeño error. La función está para debuggear.
    """

    # Cálculo de la longitud
    lon = math.atan2(y, x)

    # Cálculo de la latitud mediante iteración de Newton
    p = math.sqrt(x ** 2 + y ** 2)
    lat = math.atan2(z, p * (1 - f))
    N = a / math.sqrt(1 - e2 * math.sin(lat) ** 2)
    prev_lat = lat
    while True:
        lat = math.atan2(z + e2 * N * math.sin(lat), p)
        if abs(lat - prev_lat) < eps:
            break
        prev_lat = lat

    # Cálculo de la altitud
    alt = p / math.cos(lat) - N

    # Convertir latitud y longitud a grados
    lat = math.degrees(lat)
    lon = math.degrees(lon)

    return lat, lon, alt


def calculate_enu_from_geodetic_ecef(observer_coords: np.ndarray, target_ecef: np.ndarray) -> Tuple[
    float, float, float]:
    """Calculates ENU coordinates from geodetic coordinates of the observer and the target.

    Args:
        observer_coords (np.ndarray): Observer's geodetic coordinates [lat, lon, alt].
        target_coords (np.ndarray): Target's geodetic coordinates [lat, lon, alt].

    Returns:
        Tuple[float, float, float]: East, North, Up coordinates as a tuple.
    """
    observer_ecef = geodetic_to_ecef(*observer_coords)

    # Compute differences in ECEF coordinates
    dx = target_ecef - observer_ecef

    # Convert latitude and longitude from degrees to radians
    lat_rad = np.radians(observer_coords[0])
    lon_rad = np.radians(observer_coords[1])

    # Calculate transformation matrix elements
    cos_lat = np.cos(lat_rad)
    sin_lat = np.sin(lat_rad)
    cos_lon = np.cos(lon_rad)
    sin_lon = np.sin(lon_rad)

    # Compute ENU coordinates
    east = -sin_lon * dx[0] + cos_lon * dx[1]
    north = -sin_lat * cos_lon * dx[0] - sin_lat * sin_lon * dx[1] + cos_lat * dx[2]
    up = cos_lat * cos_lon * dx[0] + cos_lat * sin_lon * dx[1] + sin_lat * dx[2]

    return (east, north, up)