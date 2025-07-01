from ollikit import Olig_Finder
import os
from pathlib import Path

input_data ={"target_seqs":["[BHQ2]AAACAGATAGTAC",
                "AAACAGATAGTAC",
            "AAACAGATAGTAC[Cy3]"],
"target_names":["s1", "s2", "s3"],
"target_aff_low":["-35000" ,"-35000", "-35000"],
"target_aff_high":["-45000", "-45000", "-45000"],
"Hairpin_energy_thr":0,
"metric":"gibbs",
"celsius":25,
"num_oligos":1,
"aff_predictor":"Nupack",
"hairpin_predictor":"Seqfold"
}
#input_data = { "target_seqs": [ "AGAAACACGGAGUUUCGCAC", "GUGCGAUCCUCGGUGAAUCU" ], "target_names": [ "s1", "s2" ], "sequence_types": [ "RNA", "RNA" ], "target_aff_low": [ 0, 0 ], "target_aff_high": [ 1, 1 ], "Hairpin_energy_thr": 0, "metric": "fraction", "celsius": 25, "num_oligos": 1, "aff_predictor": "Nupack", "hairpin_predictor": "Seqfold", "output_format": "DNA" }
input_data = {"target_seqs":["ACGCCGCCTCAATAGATCCTACGCCGCCTC","GAGGCGGCGT","AATCCGCCTC"],"target_names":["DNAtemplate","I (product)","Q1min"],"sequence_types":["DNA","DNA","DNA"],"target_aff_low":[0,0,0.88],"target_aff_high":[0.1,0.1,0.93],"Hairpin_energy_thr":-3,"metric":"fraction","celsius":25,"num_oligos":1,"aff_predictor":"Nupack","hairpin_predictor":"Nupack","output_format":"DNA","timeout":84600}

def run_olig_finder(input_file, output_folder):
    of = Olig_Finder(input_file, output_folder)
    res = of.find_olig_robust()
    df = of.describe_oligos(res)
    of.save_to_excel(df)
    print("Olig_Finder: \n", df)

output_folder = "test_output"
os.makedirs(output_folder, exist_ok=True)  # Создаём папку, если её нет
run_olig_finder(input_data, "output")