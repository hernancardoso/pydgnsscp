import numpy as np

def compute_receiver_position_least_squares(satellites, initial_guess):
    """
    Computes the receiver's position using the Least Squares method without iterations.

    :param satellites: List of dictionaries, each containing:
                       - 'x', 'y', 'z': Satellite ECEF coordinates (meters)
                       - 'pseudorange': Measured pseudorange (meters)
    :return: Dictionary with receiver's position and clock bias:
             - 'x', 'y', 'z': Receiver ECEF coordinates (meters)
             - 'delta_t': Receiver clock bias (seconds)
    """
    c = 299792458  # Speed of light in m/s

    num_sats = len(satellites)
    if num_sats < 4:
        return None #raise ValueError("At least four satellites are required for a 3D position fix.")

    # Extract satellite positions and pseudoranges
    positions = np.array([[sat['x'], sat['y'], sat['z']] for sat in satellites])
    pseudoranges = np.array([sat['pseudorange'] for sat in satellites])

    # Choose a reference satellite (first one)
    x_ref, y_ref, z_ref = initial_guess
    p_ref = pseudoranges[0]

    # Formulate the linear system
    A = np.zeros((num_sats - 1, 4))
    b = np.zeros(num_sats - 1)

    for i in range(1, num_sats):
        x_i, y_i, z_i = positions[i]
        p_i = pseudoranges[i]

        # Compute differences
        delta_p = p_i - p_ref
        delta_x = x_i - x_ref
        delta_y = y_i - y_ref
        delta_z = z_i - z_ref

        # Approximate range differences (assuming small changes)
        A[i - 1, 0] = delta_x
        A[i - 1, 1] = delta_y
        A[i - 1, 2] = delta_z
        A[i - 1, 3] = c * (p_i - p_ref)

        b[i - 1] = 0.5 * (delta_p**2 - delta_x**2 - delta_y**2 - delta_z**2)

    # Solve the linear system using least squares
    # x = [x_receiver, y_receiver, z_receiver, c * delta_t]
    x, residuals, rank, s = np.linalg.lstsq(A, b, rcond=None)

    # Compute the receiver's position
    x_receiver = x_ref + x[0]
    y_receiver = y_ref + x[1]
    z_receiver = z_ref + x[2]
    delta_t = x[3] / (c**2)

    result = {
        'x': x_receiver,
        'y': y_receiver,
        'z': z_receiver,
        'delta_t': delta_t  # Receiver clock bias in seconds
    }
    return result

# Example usage
if __name__ == "__main__":
    # Example satellite data
    satellites = [
        {'x': 15600000.0, 'y': 7540000.0,  'z': 20140000.0, 'pseudorange': 21407000.0},
        {'x': 18760000.0, 'y': 2750000.0,  'z': 18610000.0, 'pseudorange': 21170000.0},
        {'x': 17610000.0, 'y': 14630000.0, 'z': 13480000.0, 'pseudorange': 22770000.0},
        {'x': 19170000.0, 'y': 610000.0,   'z': 18610000.0, 'pseudorange': 21380000.0},
        # Add more satellites if available
    ]

    # Compute the receiver's position
    receiver_pos = compute_receiver_position_least_squares(satellites)

    # Print the results
    print("Estimated Receiver Position and Clock Bias:")
    print(f"x: {receiver_pos['x']:.3f} meters")
    print(f"y: {receiver_pos['y']:.3f} meters")
    print(f"z: {receiver_pos['z']:.3f} meters")
    print(f"Clock Bias: {receiver_pos['delta_t']:.9f} seconds")
