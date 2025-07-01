from ollikit import Gene_Olig_Finder
import os
from pathlib import Path

def run_gene_oligfinder(input_data, output_folder):
    """
    Запускает подбор цепочки олигонуклеотидов и выводит результаты.
    """
    # Создаём экземпляр класса Gene_Olig_Finder
    finder = Gene_Olig_Finder(input_data, output_folder)
    
    # Ищем цепочку олигонуклеотидов
    print("Запуск Gene_Olig_Finder...")
    oligos = finder.find_olig_joblib()

    # Сохраняем и описываем результаты
    affinity_df, params_df, input_params_df, target_df = finder.describe_oligos(oligos)
    finder.save_to_excel(affinity_df, params_df, input_params_df, target_df)
    print("Результаты Gene_Olig_Finder:\n", affinity_df)

    print(f"Результаты сохранены в папке: {output_folder}")


input_data = {
    "target_seqs": [
        "AGAAACAC"
    ],
    "target_names": [
        "s1"
    ],
    "target_aff_low": [
        0
    ],
    "target_aff_high": [
        1
    ],
    "Hairpin_energy_thr": 0,
    "metric": "fraction",
    "celsius": 25,
    "num_oligos": 1,
    "aff_predictor": "Nupack",
    "hairpin_predictor": "Seqfold",
    "gene": "TACTGATCGATGCATCGTAGCATCGTACGATCGTAGCTGACG",
    "site_start": 18,
    "site_length": 10,
    "max_consecutive_ratio": 0.4,
    "max_overall_ratio": 0.3,
    "gen_aff_low": 0.5,
    "gen_aff_high": 1
}

# Указываем папку для сохранения результатов
output_folder = "test_output"
os.makedirs(output_folder, exist_ok=True)  # Создаём папку, если её нет

#input_data = "/biograph/prod/storage/tasks/task_3524.json"
# Запуск тестового скрипта
if __name__ == "__main__":
    # Подключаем класс Gene_Olig_Finder
    run_gene_oligfinder(input_data, output_folder)
