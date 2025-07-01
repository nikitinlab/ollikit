import logging
import numpy as np
import os
from pathlib import Path
from ollikit import Affinity_Matrix

# Настройка логирования
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('affinity_predictor_test.log'),
        logging.StreamHandler()
    ]
)

# Тестовые последовательности
test_seqs = [
    "ACTGCTAGAGATTTTCCACAT",
    "TGTGGAAAATCTCTAGCAGTT",
    "AATCGCTAGCTAGCTAGCTAG"
]

# Пример входных данных
input_data = {
    "target_seqs": test_seqs,
    "target_names": ["s1", "s2", "s3"],
    "metric": "fraction",
    "celsius": 25,
    "aff_predictor": "Nupack",
    "hairpin_predictor": "Seqfold"
}
output_folder = Path("test_output/test_affmatrix")
output_folder.mkdir(parents=True, exist_ok=True)

def test_affinity_matrix():
    aff_matr_obj = Affinity_Matrix(input_data, output_folder)
    aff_matr = aff_matr_obj.predict()  # Получаем матрицу аффинности (numpy array)
    print("Affinity_Matrix:\n", aff_matr)
    assert isinstance(aff_matr, np.ndarray)
    assert aff_matr.shape[0] == aff_matr.shape[1]  # Квадратная матрица

    # Сохраняем heatmap и Excel-файл
    aff_matr_obj.heatmaps_to_png(aff_matr)
    aff_matr_obj.matrix_to_excel(aff_matr)
    print("Heatmaps и Excel-файл успешно сохранены.")

    assert os.path.exists(f"{output_folder}/Affinity_heatmap.png")
    assert os.path.exists(f"{output_folder}/Hairpin_heatmap.png")
    assert os.path.exists(f"{output_folder}/Affinity_matrix.xlsx")

if __name__ == "__main__":
    test_affinity_matrix()