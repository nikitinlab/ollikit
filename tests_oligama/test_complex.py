from ollikit import Complex_Solver
from pathlib import Path
import os

def run_complex_solver(input_file, output_folder):
    solver = Complex_Solver(input_file, output_folder)
    data = solver.predict()
    solver.save_to_excel(data)
    solver.save_plot(data)
    print("Complex_Solver: \n", data)


input = {
            "target_seqs":  
                [
                "TGCACTATGGCACACTGGTA",
                "AGATTCGCCGTTAATCGCAA",
                "TCGAATTCCATTGTGCCATA",
                "AATAAGCGAGGCAGAGAATA",
                "AATGATCTAGCCCGCCGTTT",
                "ACCGAGAGTATTCAGCTTTA"
                ],
            "target_names":["s1", "s2", "s3", "s4", "s5", "s6"],
            "sequence_type":["DNA", "RNA", "DNA", "DNA", "RNA", "DNA"],
            "total_conc":["1e-6", "1e-6", "1e-6", "1e-6", "1e-6", "1e-6"],
            "metric":"fraction",
            "celsius":25,
            "aff_predictor":"Nupack",
            "hairpin_predictor":"Seqfold"
        }
output_folder = Path("test_output/test_complex")
output_folder.mkdir(parents=True, exist_ok=True)

run_complex_solver(input, "output")
