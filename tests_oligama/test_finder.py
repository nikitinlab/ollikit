from ollikit import Olig_Finder
import os
from pathlib import Path
import tempfile
import pandas as pd
import pytest

def run_olig_finder(input_file, output_folder):
    of = Olig_Finder(input_file, output_folder)
    res = of.find_olig_robust()
    df = of.describe_oligos(res)
    of.save_to_excel(df)
    print("Olig_Finder: \n", df)

def test_olig_finder():
    input_data = { "target_seqs": [ "AGAAACACGGAGUUUCGCAC", "GUGCGAUCCUCGGUGAAUCU" ], "target_names": [ "s1", "s2" ], "sequence_types": [ "RNA", "RNA" ], "target_aff_low": [ 0, 0 ], "target_aff_high": [ 1, 1 ], "Hairpin_energy_thr": 0, "metric": "fraction", "celsius": 25, "num_oligos": 1, "aff_predictor": "Nupack", "hairpin_predictor": "Seqfold", "output_format": "RNA" }

    with tempfile.TemporaryDirectory() as tmpdir:
        finder = Olig_Finder(input_data, tmpdir)
        res = finder.find_olig_robust()
        df = finder.describe_oligos(res)
        assert isinstance(df, pd.DataFrame)
        out_file = os.path.join(tmpdir, "olig_finder_test.xlsx")
        finder.save_to_excel(df, filename="olig_finder_test.xlsx")
        assert os.path.exists(out_file)
        print("Olig_Finder DataFrame:\n", df)

# output_folder = "test_output"
# os.makedirs(output_folder, exist_ok=True)  # Создаём папку, если её нет
# run_olig_finder(input_data, "output")