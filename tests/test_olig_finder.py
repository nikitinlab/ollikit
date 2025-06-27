import tempfile
import pandas as pd
from ollikit import Olig_Finder

def test_olig_finder_basic():
    input_data = {
        "target_seqs": ["ACGCCGCCTCAATAGATCCTACGCCGCCTC", "GAGGCGGCGT", "AATCCGCCTC"],
        "target_names": ["DNAtemplate", "I (product)", "Q1min"],
        "sequence_types": ["DNA", "DNA", "DNA"],
        "target_aff_low": [0, 0, 0.88],
        "target_aff_high": [0.1, 0.1, 0.93],
        "Hairpin_energy_thr": -3,
        "metric": "fraction",
        "celsius": 25,
        "num_oligos": 1,
        "aff_predictor": "Nupack",
        "hairpin_predictor": "Nupack",
        "output_format": "DNA",
        "timeout": 10
    }
    with tempfile.TemporaryDirectory() as tmpdir:
        finder = Olig_Finder(input_data, tmpdir)
        res = finder.find_olig_robust()
        df = finder.describe_oligos(res, toxls=False)
        assert isinstance(df, pd.DataFrame) 