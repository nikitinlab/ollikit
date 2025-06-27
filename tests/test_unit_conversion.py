import numpy as np
import pytest
from ollikit import triangle_to_matr, detect_unit, convert

def test_triangle_to_matr():
    vec = [1, 2, 3]
    matr = triangle_to_matr(vec)
    assert matr.shape == (3, 3)
    assert np.allclose(matr, matr.T)
    assert np.allclose(np.triu(matr, 1)[np.triu_indices(3, 1)], vec)

def test_detect_unit():
    assert detect_unit([-5, -2, -1]) == "energy"
    assert detect_unit([2, 3, 4]) == "Ka"
    assert detect_unit([1e-5, 0.01, 0.05]) == "Kd"
    assert detect_unit([0.2, 0.5, 0.8]) == "fraction"

def test_convert_fraction_to_ka():
    arr = np.array([0.1, 0.5])
    ka = convert(arr, to_unit="ka", from_unit="fraction")
    assert ka.shape == arr.shape
    assert np.all(ka > 0)

def test_convert_ka_to_fraction():
    arr = np.array([1e4, 1e5])
    fraction = convert(arr, to_unit="fraction", from_unit="ka")
    assert fraction.shape == arr.shape
    assert np.all((fraction >= 0) & (fraction <= 1)) 