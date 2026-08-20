import numpy as np
import tempfile
import os

from class_cam_generation import read_xy_from_csv, CamGeneration


def test_read_xy_from_csv(tmp_path):
    p = tmp_path / "data.csv"
    p.write_text("1,2\n3,4\n")
    x, y = read_xy_from_csv(str(p))
    assert np.allclose(x, np.array([1.0, 3.0]))
    assert np.allclose(y, np.array([2.0, 4.0]))


def test_remove_duplicates_and_to_polar():
    # Create a minimal CamGeneration instance
    gear_ratios = np.array([1.0, 1.0])
    input_angles = np.array([0.0, np.pi])
    cam = CamGeneration(gear_ratios, input_angles, scaling=1.0, sit_angle=0.0, offset_angle=0.0)

    # Test remove_duplicates
    x = np.array([1.0, 2.0, 2.0, 3.0])
    y = np.array([10.0, 20.0, 20.0, 30.0])
    z = np.array([100.0, 200.0, 200.0, 300.0])
    x2, y2, z2 = cam.remove_duplicates(x, y, z)
    # remove_duplicates removes values that appear more than once entirely
    # so only unique, non-repeated values remain
    assert np.allclose(x2, np.array([1.0, 3.0]))
    assert np.allclose(y2, np.array([10.0, 30.0]))
    assert np.allclose(z2, np.array([100.0, 300.0]))

    # Test to_polar
    points = np.array([[1.0, 0.0], [0.0, 1.0], [-1.0, 0.0]])
    radii, angles = cam.to_polar(points)
    assert np.allclose(radii, np.array([1.0, 1.0, 1.0]))
    # angles should be 0, pi/2, pi (or equivalent in [0, 2*pi])
    assert np.allclose(angles, np.array([0.0, np.pi/2, np.pi]))
