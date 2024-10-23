import georinex as gr
import numpy as np
from datetime import datetime
from math import atan2, sqrt, sin, cos, degrees, radians

from dgnsscp.coordinate_transformations import ecef_to_enu


def llh2ecef(lat, lon, h):
    """Convert geodetic coordinates to ECEF."""
    # WGS84 ellipsoid constants
    a = 6378137.0  # Semi-major axis
    f = 1 / 298.257223563  # Flattening
    e2 = f * (2 - f)  # Square of eccentricity

    lat_rad = radians(lat)
    lon_rad = radians(lon)

    N = a / sqrt(1 - e2 * sin(lat_rad) ** 2)

    x = (N + h) * cos(lat_rad) * cos(lon_rad)
    y = (N + h) * cos(lat_rad) * sin(lon_rad)
    z = (N * (1 - e2) + h) * sin(lat_rad)

    return x, y, z


def calculate_az_el(sat_pos_ecef, receiver_pos_llh):
    """Calculate azimuth and elevation from receiver to satellite."""
    # Satellite position in ECEF
    x, y, z = sat_pos_ecef

    # Receiver position in LLH
    lat0, lon0, h0 = receiver_pos_llh

    # Convert satellite position to ENU coordinates
    e, n, u = ecef_to_enu(x, y, z, lat0, lon0, h0)

    # Calculate azimuth and elevation
    azimuth = atan2(e, n)
    elevation = atan2(u, sqrt(e ** 2 + n ** 2))

    # Convert from radians to degrees
    azimuth_deg = (degrees(azimuth) + 360) % 360
    elevation_deg = degrees(elevation)

    return azimuth_deg, elevation_deg


def compute_satellite_position(eph, obs_time):
    """Compute satellite ECEF position from ephemeris at observation time."""
    # Constants
    mu = 3.986005e14  # Earth's universal gravitational parameter (m^3/s^2)
    Omega_e_dot = 7.2921151467e-5  # Earth's rotation rate (rad/s)

    # Time from ephemeris reference epoch
    toe = eph['Toe'].values.item()
    t = (obs_time - datetime(1980, 1, 6)).total_seconds()
    tk = t - toe

    # Correct for end-of-week crossover
    if tk > 302400:
        tk -= 604800
    elif tk < -302400:
        tk += 604800

    # Extract ephemeris parameters
    sqrtA = eph['sqrtA'].values.item()
    A = sqrtA ** 2
    e = eph['Eccentricity'].values.item()
    i0 = eph['Io'].values.item()
    Omega0 = eph['Omega0'].values.item()
    omega = eph['omega'].values.item()
    M0 = eph['M0'].values.item()
    Delta_n = eph['DeltaN'].values.item()
    IDOT = eph['IDOT'].values.item()
    Cuc = eph['Cuc'].values.item()
    Cus = eph['Cus'].values.item()
    Crc = eph['Crc'].values.item()
    Crs = eph['Crs'].values.item()
    Cic = eph['Cic'].values.item()
    Cis = eph['Cis'].values.item()
    Omega_dot = eph['OmegaDot'].values.item()

    # Corrected mean motion
    n0 = sqrt(mu / A ** 3)
    n = n0 + Delta_n

    # Mean anomaly
    M = M0 + n * tk

    # Solve Kepler's Equation for Eccentric Anomaly E
    E = M
    for _ in range(10):
        E_old = E
        E = M + e * sin(E)
        if abs(E - E_old) < 1e-12:
            break

    # True anomaly
    sin_v = (sqrt(1 - e ** 2) * sin(E)) / (1 - e * cos(E))
    cos_v = (cos(E) - e) / (1 - e * cos(E))
    v = atan2(sin_v, cos_v)

    # Argument of latitude
    phi = v + omega

    # Second harmonic perturbations
    du = Cuc * cos(2 * phi) + Cus * sin(2 * phi)
    dr = Crc * cos(2 * phi) + Crs * sin(2 * phi)
    di = Cic * cos(2 * phi) + Cis * sin(2 * phi)

    # Corrected parameters
    u = phi + du
    r = A * (1 - e * cos(E)) + dr
    i = i0 + di + IDOT * tk

    # Positions in orbital plane
    x_prime = r * cos(u)
    y_prime = r * sin(u)

    # Corrected longitude of ascending node
    Omega = Omega0 + (Omega_dot - Omega_e_dot) * tk - Omega_e_dot * toe

    # ECEF coordinates
    x = x_prime * cos(Omega) - y_prime * cos(i) * sin(Omega)
    y = x_prime * sin(Omega) + y_prime * cos(i) * cos(Omega)
    z = y_prime * sin(i)

    return np.array([x, y, z])


def calculate_az_el_dict(receiver_pos_llh, obs_time, navfile):
    # Read RINEX navigation file
    nav = gr.load(navfile)
    az_el_dict = {}
    # Get list of satellites
    svs = nav.sv.values

    # print(f"Calculating Azimuth and Elevation for {obs_time} UTC")
    # print(f"Receiver Position: Lat {receiver_pos_llh[0]}, Lon {receiver_pos_llh[1]}, Height {receiver_pos_llh[2]} m")
    # print("Satellite      Azimuth (deg)  Elevation (deg)")

    for sv in svs:
        try:
            # Get satellite ephemeris at the closest time
            sat_data = nav.sel(sv=sv).dropna(dim='time', how='all')
            times = sat_data.time.values

            # Convert obs_time to numpy.datetime64
            obs_time_np = np.datetime64(obs_time)

            # Find the closest time index
            time_deltas = np.abs(times - obs_time_np)
            idx = np.argmin(time_deltas)

            # Get ephemeris parameters
            eph = sat_data.isel(time=idx)

            # Compute satellite position in ECEF coordinates
            sat_pos_ecef = compute_satellite_position(eph, obs_time)

            # Calculate azimuth and elevation
            azimuth_deg, elevation_deg = calculate_az_el(sat_pos_ecef, receiver_pos_llh)
            sat_prn = int(sv[1:])
            az_el_dict[sat_prn] = (azimuth_deg, elevation_deg)

            # print(f"{sv}        {azimuth_deg:14.2f}  {elevation_deg:14.2f}")
        except Exception as e:
            print(f"{sv}        Calculation Error: {e}")
    return az_el_dict


if __name__ == '__main__':
    # Receiver fixed position (latitude, longitude, height in meters)
    receiver_pos_llh = (0, 0, 0)  # Replace with your coordinates

    # Observation time (UTC)
    obs_time = datetime(2024, 10, 1, 5, 0, 0)  # Replace with your desired time
    navfile = 'brdc2680.24n'  # Replace with your file name

    calculate_az_el_dict(receiver_pos_llh, obs_time, navfile)
