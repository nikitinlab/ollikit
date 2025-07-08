import numpy as np
import pytest
from ollikit.oligama.predictors.affinity_predictors import Nupack_Affinity_Predictor

def test_predict_two_strings():
    seq1 = "ACTGCTAGAGATTTTCCACAT"
    seq2 = "TGTGGAAAATCTCTAGCAGTT"
    predictor = Nupack_Affinity_Predictor(celsius=25, material='dna')
    result = predictor.predict(seq1, seq2, units='fraction')
    print("Результат аффинности для двух строк:", result)
    assert isinstance(result, np.ndarray)
    assert result.shape[0] == 1 or result.shape == ()  # допускаем скаляр или массив из одного элемента

def test_predict_arrays():
    seqs1 = ["ACTGCTAGAGATTTTCCACAT", "AATCGCTAGCTAGCTAGCTAG"]
    seqs2 = ["TGTGGAAAATCTCTAGCAGTT", "ACTGCTAGAGATTTTCCACAT"]
    predictor = Nupack_Affinity_Predictor(celsius=25, material='dna')
    result_arr = predictor.predict(seqs1, seqs2, units='fraction')
    print("Результаты аффинности для массивов строк:", result_arr)
    assert isinstance(result_arr, np.ndarray)
    assert result_arr.shape[0] == 2 