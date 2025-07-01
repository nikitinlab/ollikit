from ollikit import Primary_Concentration_Solver
from pathlib import Path
import os

def run_concentration_solver(input_file, output_folder):
    solver = Primary_Concentration_Solver(input_file, output_folder)
    df, inp_df = solver.predict()
    solver.save_plot(df)
    solver.save_to_excel(df, inp_df)
    print("Concentration_Solver: \n", df)

input = {"target_seqs":["AGAAACACGGAGTTTCGCAC","GTGCGATCCTCGGTGAATCT","AGAATCACAGAGGATCGCAG"],"target_names":["s1","s2","s3"],"sequence_types":["DNA","DNA","DNA"],"total_conc":["1e-6","1e-6","1e-6"],"input_seq_name":"s1","output_seq_name":"s2","input_seq_bounds":["1e-9","1e-6"],"n_points":20,"aff_predictor":"Oligama","hairpin_predictor":"Seqfold","metric":"fraction","celsius":25}

output_folder = Path("test_output/test_conc")
os.mkdir(parents=True, exist_ok=True)

run_concentration_solver(input, output_folder)
