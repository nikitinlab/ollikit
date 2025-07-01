import tempfile
import pandas as pd
from ollikit import Primary_Concentration_Solver
import os

def test_primary_concentration_solver_basic():
    input_data = {
        "target_seqs": ["ACGCCGCCTCAATAGATCCTACGCCGCCTC", "GAGGCGGCGT", "AATCCGCCTC"],
        "target_names": ["DNAtemplate", "I (product)", "Q1min"],
        "sequence_types": ["DNA", "DNA", "DNA"],
        "total_conc": ["1e-6", "1e-6", "1e-6"],
        "input_seq_name": "DNAtemplate",
        "output_seq_name": "I (product)",
        "input_seq_bounds": [1e-8, 1e-5],
        "n_points": 5,
        "metric": "fraction",
        "celsius": 25,
        "aff_predictor": "Nupack",
        "hairpin_predictor": "Nupack",
        "output_format": "DNA"
    }
    with tempfile.TemporaryDirectory() as tmpdir:
        solver = Primary_Concentration_Solver(input_data, tmpdir)
        df, inp_df = solver.predict()
        assert isinstance(df, pd.DataFrame)
        assert not df.empty
        solver.save_to_excel(df, inp_df, filename="primary_conc_solver_test.xlsx")
        assert os.path.exists(os.path.join(tmpdir, "primary_conc_solver_test.xlsx"))