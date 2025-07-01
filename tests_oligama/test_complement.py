from ollikit import Complement
from pathlib import Path
import os

def run_complement(input_file, output_folder):
    comp = Complement(input_file, output_folder)
    rev_arr, rev_compl_arr, dir_compl_arr = comp.predict()
    comp.save_to_excel(rev_arr, rev_compl_arr, dir_compl_arr)
    print("Complement: \n", rev_arr, rev_compl_arr, dir_compl_arr)

input = {"target_seqs":["AGAAACACGGAGUUUCGCAC"],"target_names":["q"],"sequence_types":["RNA"]}

output_folder = Path("test_output/test_complement")
output_folder.mkdir(parents=True, exist_ok=True)
run_complement(input, output_folder)
