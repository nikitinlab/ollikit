import numpy as np
import pytest
from ollikit import Nupack_Affinity_Predictor, Nupack_Hairpin_Predictor

def test_nupack_affinity_predictor_fraction():
    predictor = Nupack_Affinity_Predictor(celsius=25)
    seq1 = ["ATGC"]
    seq2 = ["GCAT"]
    res = predictor.predict(seq1, seq2, units="fraction")
    assert isinstance(res, np.ndarray)
    assert res.shape == (1,)
    assert np.all(res >= 0)

def test_nupack_affinity_predictor_kd():
    predictor = Nupack_Affinity_Predictor(celsius=25)
    seq1 = ["ATGC"]
    seq2 = ["GCAT"]
    res = predictor.predict(seq1, seq2, units="Kd")
    assert isinstance(res, np.ndarray)
    assert res.shape == (1,)
    assert np.all(res > 0)

def test_nupack_hairpin_predictor():
    predictor = Nupack_Hairpin_Predictor(celsius=25)
    seqs = ["ATGCATGC", "GGGGCCCC"]
    res = predictor.predict(seqs)
    assert isinstance(res, np.ndarray)
    assert res.shape == (2,) 