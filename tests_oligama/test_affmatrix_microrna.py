import os
import pandas as pd
from ollikit import Affinity_Matrix  # Импортируем нужный класс
from pathlib import Path

# данные для тестирования MicroRNA_Finder
input_data = {
                "target_seqs":["UGAGGUAGUAGGUUGUGUGGUU",
                "CUGUACAACCUUCUAGCUUUCC"],
                "target_names":["s", "i"],
                "metric":"fraction",
                "celsius":25,
                "aff_predictor":"Oligama",
                "hairpin_predictor":"Seqfold"
                }

def update_target_fields(input_data, microrna_dict):
    """
    Обновляет target_seqs и target_names в input_data на основе microrna_dict.
    
    Args:
        input_data (dict): Исходные данные, содержащие target_seqs и target_names.
        microrna_dict (dict): Словарь microRNA, где ключи — это имена, а значения — последовательности.
    
    Returns:
        dict: Обновленные данные input_data с новыми target_seqs и target_names.
    """
    input_data["target_seqs"] = list(microrna_dict.values())  # Берем значения из словаря
    input_data["target_names"] = list(microrna_dict.keys())   # Берем ключи из словаря
    return input_data

def load_microrna_list():
    """Загружает microRNA из файла."""
    base_path = os.path.dirname(__file__)
    file_path = os.path.join(base_path, "oligama/helpers/microrna/mature.txt")
    A = pd.read_csv(file_path, sep='\t', header=None).to_numpy()
    miRNA = {}
    for i in range(len(A) // 2):
        if 'Homo sapiens' in str(A[i]):
            key = str(A[i]).split(' ')[0][3:]
            seq = str(A[i + 1])[2:-2]
            miRNA[key] = seq
    return miRNA

def run_affinity_matrix(input_data, output_folder):
    aff_matr = Affinity_Matrix(input_data, output_folder)
    aff_matr.predict()
    print("Affinity_Matrix: \n", aff_matr.predict())

# Загружаем microRNA
#microrna_dict = load_microrna_list()
# Обновляем target_seqs и target_names в input_data
#input_data = update_target_fields(input_data, microrna_dict)

# Печатаем обновленный input_data
print("Updated input_data:")
print(input_data)
output_folder = Path("test_output/test_affmatrix_microrna")
output_folder.mkdir(parents=True, exist_ok=True)

run_affinity_matrix(input_data, output_folder)