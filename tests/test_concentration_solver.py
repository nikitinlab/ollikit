import numpy as np
from ollikit import fpi_step, find_equilibrium_conc

def test_fpi_step():
    total_conc = np.array([1.0, 2.0], dtype=np.float64)
    aff_matr_ka = np.array([[0, 1], [1, 0]], dtype=np.float64)
    x = np.array([0.5, 1.5], dtype=np.float64)
    res = fpi_step(x, aff_matr_ka, total_conc)
    assert res.shape == x.shape
    assert np.all(res >= 0)

def test_find_equilibrium_conc():
    total_conc = np.array([1.0, 2.0], dtype=np.float64)
    aff_matr_ka = np.array([[0, 1], [1, 0]], dtype=np.float64)
    res = find_equilibrium_conc(total_conc, aff_matr_ka)
    assert res.shape == total_conc.shape
    assert np.all(res >= 0) 