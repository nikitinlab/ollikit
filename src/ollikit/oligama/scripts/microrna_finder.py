import pandas as pd
import numpy as np
from ..utils import df_to_excel
from ..data_loaders import MicroRNA_Dataloader
import os


class MicroRNA_Finder(MicroRNA_Dataloader):
    def __init__(self, input_data, output_folder):
        super().__init__(input_data, output_folder)
        # Дополнительные параметры
        self.num_micrornas = self.data.get("num_micrornas", 5)  # Количество требуемых microRNA
        self.microrna_dict = self.load_microrna_list()  # Загрузка списка microRNA

    def load_microrna_list(self):
        """Загружает microRNA из файла."""
        base_path = os.path.dirname(__file__)
        file_path = os.path.join(base_path, "../helpers/microrna/mature.txt")
        A = pd.read_csv(file_path, sep='\t', header=None).to_numpy()
        miRNA = {}
        for i in range(len(A) // 2):
            if 'Homo sapiens' in str(A[i]):
                key = str(A[i]).split(' ')[0][3:]
                seq = str(A[i + 1])[2:-2]
                miRNA[key] = seq
        return miRNA

    def batch_affinity_predict(self, seq1_arr, seq2_arr):
        """Вычисляет афинность с использованием предиктора."""
        master_seq_arr = self.arr_product(seq1_arr, seq2_arr)
        aff_arr = self.aff_predictor.predict(master_seq_arr[:, 0], master_seq_arr[:, 1], units=self.metric)
        aff_arr = aff_arr.reshape((len(seq1_arr), len(seq2_arr)))
        return aff_arr

    def arr_product(self, arr1, arr2):
        """Создает все возможные пары последовательностей для расчета афинности."""
        len1, len2 = len(arr1), len(arr2)
        arr1 = np.array([arr1] * len2).reshape((len2, len1)).T
        arr2 = np.array([arr2] * len1).reshape((len1, len2))
        return np.concatenate((arr1[:, :, None], arr2[:, :, None]), axis=2).reshape(-1, 2)

    def find_best_microrna(self, input_oligs):
        """Находит наиболее комплементарные microRNA для всех входных олигов."""
        results = []
        microrna_names = list(self.microrna_dict.keys())
        microrna_seqs = list(self.microrna_dict.values())

        # Рассчитываем афинности для всех пар input_oligs и microRNA
        affinity_matrix = self.batch_affinity_predict(input_oligs, microrna_seqs)

        # Определяем порядок сортировки в зависимости от self.metric
        reverse_sort = not (self.metric in ['Kd', 'gibbs'])  # True - по убыванию, False - по возрастанию

        # Для каждого входного олига выбираем top-N microRNA с наибольшей афинностью
        for i, olig in enumerate(input_oligs):
            affinities = affinity_matrix[i]
            # Сортировка с учётом направления (прямой или обратный порядок)
            top_indices = np.argsort(affinities)[:self.num_micrornas] if not reverse_sort else \
                        np.argsort(affinities)[-self.num_micrornas:][::-1]
            
            for idx in top_indices:
                results.append({
                    "Input Olig": olig,
                    "miRNA Name": microrna_names[idx],
                    "miRNA Sequence": microrna_seqs[idx],
                    "Affinity": affinities[idx]
                })

        return results

    def run(self):
        """Запускает процесс поиска лучших microRNA и возвращает результат как DataFrame."""
        results = self.find_best_microrna(self.target_seqs)
        results_df = pd.DataFrame(results)
        return results_df

    def save_to_excel(self, results_df):
        output_path = self.output_folder / "MicroRNA_Finder_Output.xlsx"
        df_to_excel([results_df], ["Best microRNA"], output_path)
