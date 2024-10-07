import numpy as np

def retrieve_ephemeris(msg):
    """
    Decode RTCM message to retrieve ephemeris data for GPS satellites.

    Parameters:
    rtcm_message (bytes): The RTCM message in bytes format.

    Returns:
    dict: A dictionary containing ephemeris data.
    """
    ephemeris = {}

    ephemeris['A'] = (msg.DF092 ** 2)  # Semi-major axis in meters (Square root of Semi-major Axis)
    ephemeris['e'] = msg.DF090  # Eccentricity
    ephemeris['i'] = msg.DF097  # Inclination in radians
    ephemeris['Ω'] = msg.DF095  # Longitude of ascending node in radians
    ephemeris['ω'] = msg.DF099  # Argument of perigee in radians
    ephemeris['M'] = msg.DF088  # Mean anomaly in radians
    ephemeris['toe'] = msg.DF093  # Time of ephemeris in seconds

    # Convert angles from semicircles to radians
    ephemeris['i'] = np.radians(ephemeris['i'])
    ephemeris['Ω'] = np.radians(ephemeris['Ω'])
    ephemeris['ω'] = np.radians(ephemeris['ω'])
    ephemeris['M'] = np.radians(ephemeris['M'])

    return ephemeris




def calculate_ecef_satellite_position(ephemeris) -> np.ndarray:
    A = ephemeris['A']
    e = ephemeris['e']
    i = ephemeris['i']
    Ω = ephemeris['Ω']
    ω = ephemeris['ω']
    M = ephemeris['M']
    toe = ephemeris['toe']
    t = toe

    # Constants
    mu = 3.986005e14  # Earth's universal gravitational parameter, m^3/s^2

    # Compute the mean motion (rad/s)
    n0 = np.sqrt(mu / A ** 3)

    # Time from ephemeris reference epoch
    tk = t - toe

    # Mean anomaly
    Mtk = M + n0 * tk

    # Solve Kepler's equation for Eccentric Anomaly (E)
    E = Mtk
    for _ in range(10):
        E = Mtk + e * np.sin(E)

    # True anomaly (ν)
    ν = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2), np.sqrt(1 - e) * np.cos(E / 2))

    # Argument of Latitude (u)
    u = ω + ν

    # Radius (r)
    r = A * (1 - e * np.cos(E))

    # Positions in orbital plane
    x_orbital = r * np.cos(u)
    y_orbital = r * np.sin(u)

    # Corrected longitude of ascending node (Ωk)
    Ωk = Ω + (0 * tk) - (7.2921151467e-5 * toe)  # Earth rotation rate correction

    # Positions in ECEF
    x = x_orbital * np.cos(Ωk) - y_orbital * np.cos(i) * np.sin(Ωk)
    y = x_orbital * np.sin(Ωk) + y_orbital * np.cos(i) * np.cos(Ωk)
    z = y_orbital * np.sin(i)

    return np.array([x, y, z])


def ecef_to_enu(receiver_lat, receiver_lon, receiver_alt, sat_ecef):
    # Convert receiver position to ECEF coordinates
    a = 6378137.0  # Earth semimajor axis (WGS84)
    f = 1 / 298.257223563  # Earth flattening (WGS84)
    e2 = f * (2 - f)  # Square of eccentricity

    # Compute N, the prime vertical radius of curvature
    N = a / np.sqrt(1 - e2 * np.sin(np.radians(receiver_lat)) ** 2)

    # Compute receiver ECEF coordinates
    receiver_ecef = np.array([
        (N + receiver_alt) * np.cos(np.radians(receiver_lat)) * np.cos(np.radians(receiver_lon)),
        (N + receiver_alt) * np.cos(np.radians(receiver_lat)) * np.sin(np.radians(receiver_lon)),
        (N * (1 - e2) + receiver_alt) * np.sin(np.radians(receiver_lat))
    ])

    # ENU conversion matrix
    cos_lat = np.cos(np.radians(receiver_lat))
    sin_lat = np.sin(np.radians(receiver_lat))
    cos_lon = np.cos(np.radians(receiver_lon))
    sin_lon = np.sin(np.radians(receiver_lon))

    R = np.array([
        [-sin_lon, cos_lon, 0],
        [-cos_lon * sin_lat, -sin_lon * sin_lat, cos_lat],
        [cos_lon * cos_lat, sin_lon * cos_lat, sin_lat]
    ])

    # Compute the vector from receiver to satellite in ECEF
    diff = sat_ecef - receiver_ecef

    # Convert to ENU coordinates
    enu = R @ diff

    return enu


def calculate_azimuth_elevation(enu):
    east, north, up = enu
    azimuth = np.degrees(np.arctan2(east, north)) % 360
    elevation = np.degrees(np.arctan2(up, np.sqrt(east ** 2 + north ** 2)))
    return azimuth, elevation
