import oligama
import os
import pandas as pd

# данные для тестирования MicroRNA_Finder
input_data = {
    "target_seqs": ["ATGCGTACG", "CGTACGTAG"], 
    "target_names":["s1", "s2"],
    "num_micrornas": 3,  # Количество лучших microRNA для подбора
    "aff_predictor": "Nupack", 
    "hairpin_predictor": "Seqfold", 
    
    "metric":"Kd",
    "celsius":25,

}

def run_microrna_finder(ol, input_file, output_folder):
    mf = ol.MicroRNA_Finder(input_file, output_folder)
    results_df = mf.run()
    
    # Отображение результатов
    print("MicroRNA_Finder Results:")
    print(results_df)

# Запуск с указанными данными
output_folder = "test_output"
os.makedirs(output_folder, exist_ok=True)  # Создаём папку, если её нет

run_microrna_finder(oligama, input_data, output_folder)
