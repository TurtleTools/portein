import numpy as np
import numba as nb


# From https://github.com/numba/numba/issues/1269#issuecomment-782189134
@nb.njit
def np_apply_along_axis(func1d, axis, arr):
    assert arr.ndim == 2
    assert axis in [0, 1]
    if axis == 0:
        result = np.empty(arr.shape[1], dtype=arr.dtype)
        for i in range(len(result)):
            result[i] = func1d(arr[:, i])
    else:
        result = np.empty(arr.shape[0], dtype=arr.dtype)
        for i in range(len(result)):
            result[i] = func1d(arr[i, :])
    return result


# Adapted from: https://stackoverflow.com/questions/2964248/optimal-rotation-of-3d-model-for-2d-projection
# and
# https://math.stackexchange.com/questions/542801/rotate-3d-coordinate-system-such-that-z-axis-is-parallel-to-a-given-vector
@nb.njit
def find_best_projection(points):
    # Recenter
    points -= np_apply_along_axis(np.mean, 0, points)
    # Calculate moment of inertia tensor
    i_xy = -np.sum(points[:, 0] * points[:, 1])
    i_yz = -np.sum(points[:, 1] * points[:, 2])
    i_xz = -np.sum(points[:, 0] * points[:, 2])
    inertia = np.array(
        [
            (np.sum(points[:, 1] ** 2 + points[:, 2] ** 2), i_xy, i_xz),
            (i_xy, np.sum(points[:, 0] ** 2 + points[:, 2] ** 2), i_yz),
            (i_xz, i_yz, np.sum(points[:, 0] ** 2 + points[:, 1] ** 2)),
        ]
    )
    # Diagonalize and find new Z axis
    values, vectors = np.linalg.eig(inertia)
    new_axis = vectors[np.argmax(values)]

    # Rotate coordinate system to new Z axis
    i, j, k = np.eye(3)
    if np.sum((new_axis - k) ** 2) < 1e-8:
        return points
    theta = np.arccos(np.dot(new_axis, k))
    b = np.cross(k, new_axis)
    b_hat = b / np.sqrt(np.sum(b ** 2))
    q0 = np.cos(theta / 2)
    q1, q2, q3 = np.sin(theta / 2) * b_hat
    quaternion = np.array(
        [
            (
                q0 ** 2 + q1 ** 2 - q2 ** 2 - q3 ** 2,
                2 * (q1 * q2 - q0 * q3),
                2 * (q1 * q3 + q0 * q2),
            ),
            (
                2 * (q2 * q1 + q0 * q3),
                q0 ** 2 - q1 ** 2 + q2 ** 2 - q3 ** 2,
                2 * (q2 * q3 - q0 * q1),
            ),
            (
                2 * (q3 * q1 - q0 * q2),
                2 * (q3 * q2 + q0 * q1),
                q0 ** 2 - q1 ** 2 - q2 ** 2 + q3 ** 2,
            ),
        ]
    )
    uvw = np.zeros((3, 3))
    uvw[0] = np.dot(quaternion, i)
    uvw[1] = np.dot(quaternion, j)
    uvw[2] = np.dot(quaternion, k)
    return np.dot(points, uvw)


# Adapted from: https://stackoverflow.com/questions/32892932/create-the-oriented-bounding-box-obb-with-python-and-numpy
@nb.njit
def rotate_to_maximize_bb_height(points):
    value, vectors = np.linalg.eig(np.cov(points, y=None, rowvar=False, bias=True))
    ar = np.dot(points, vectors)
    mina = np_apply_along_axis(np.min, 0, ar)
    maxa = np_apply_along_axis(np.max, 0, ar)
    diff = (maxa - mina) * 0.5
    center = mina + diff
    corners = np.array(
        [
            (center[0] - diff[0], center[1] - diff[1]),
            (center[0] + diff[0], center[1] - diff[1]),
            (center[0] + diff[0], center[1] + diff[1]),
            (center[0] - diff[0], center[1] + diff[1]),
            (center[0] - diff[0], center[1] - diff[1]),
        ]
    )
    corners = np.dot(corners, vectors.T)
    side_1 = corners[0] - corners[1]
    side_2 = corners[1] - corners[2]
    if np.sum(side_1 ** 2) > np.sum(side_2 ** 2):
        theta = np.arctan2(side_1[0], side_1[1])
    else:
        theta = np.arctan2(side_2[0], side_2[1])

    c, s = np.cos(theta), np.sin(theta)
    rotation = np.array(((c, -s), (s, c)))
    return np.dot(points, rotation.T)
