import numpy as np


class LatLonAltCoord:
    def __init__(self, lat: float, lon: float, alt: float) -> None:
        self.lat = lat
        self.lon = lon
        self.alt = alt

    def to_numpy(self) -> np.ndarray:
        return np.array([self.lat, self.lon, self.alt])


class XYZCoord:
    def __init__(self, x: float, y: float, z: float) -> None:
        self.x = x
        self.y = y
        self.z = z

    def to_numpy(self) -> np.ndarray:
        return np.array([self.x, self.y, self.z])


class AzElCoord:
    def __init__(self, azimuth: float, elevation: float, distance: float = 0) -> None:
        self.azimuth = azimuth
        self.elevation = elevation
        self.distance = distance


