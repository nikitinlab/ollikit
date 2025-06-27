import os
import tempfile
import pandas as pd
from ollikit import Complex_Solver

def test_complex_solver_basic():
    # Минимальный пример входных данных
    input_data = {
        "target_seqs": [
            "TGCACTATGGCACACTGGTA",
            "AGATTCGCCGTTAATCGCAA",
            "TCGAATTCCATTGTGCCATA"
        ],
        "target_names": ["s1", "s2", "s3"],
        "total_conc": ["1e-6", "1e-6", "1e-6"],
        "metric": "fraction",
        "celsius": 25,
        "aff_predictor": "Nupack",
        "hairpin_predictor": "Seqfold"
    }
    # Временная папка для вывода
    with tempfile.TemporaryDirectory() as tmpdir:
        solver = Complex_Solver(input_data, tmpdir)
        result = solver.predict(toxls=False)
        assert isinstance(result, pd.DataFrame)
        assert not result.empty 