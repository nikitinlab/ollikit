import time
import numpy as np
import pandas as pd
from oligama.data_loaders import Olig_Finder_Dataloader
from oligama.utils import df_to_excel, multiple_crossover
from oligama.exceptions import OligamaWarning
import matplotlib.pyplot as plt
from joblib import Parallel, delayed
from multiprocessing import Manager
import math
    
class Gene_Olig_Finder(Olig_Finder_Dataloader):
    
    def __init__(self, input_data, output_folder):
        
        super().__init__(input_data, output_folder)
        
        self.batch_size = 360
        # self.max_iter = 40 защита от зацикливания заменена на time_limit
        self.time_limit = 600
        self.gene = self.data.get("gene", "")

        self.site_start = self.data.get("site_start", 0)
        self.site_length = self.data.get("site_length", 10)
        # Максимальное соотношение количества  комплементарных пар к длине сайта
        self.max_consecutive_ratio = self.data.get("max_consecutive_ratio", 0.3)
        # Максимальное соотношение количество подряд идущих комплементарных пар к длине сайта
        self.max_overall_ratio = self.data.get("max_overall_ratio", 0.3)




        self.gen_aff_low = self.data.get("gen_aff_low", 0)
        self.gen_aff_high = self.data.get("gen_aff_high", 0)
        
        # Рассчитываем ограничения на основе длины сайта
        # Максимальное количество подряд идущих комплементарных пар
        self.max_consecutive_letters = max(1, math.floor(self.site_length * self.max_consecutive_ratio))
        self.max_overall_letters = max(1, math.floor(self.site_length * self.max_overall_ratio))

        # Проверка порогов энергии шпилек
        target_hairpin_en = np.array(self.hairpin_predictor.predict(self.target_seqs))
        if len(target_hairpin_en) > 0 and min(target_hairpin_en) < self.hairpin_en_thr:
            mes = f"""Hairpin energy threshold ({self.hairpin_en_thr:.2f} kJ/mol) is greater than hairpin energy for 
            {self.target_names[np.argmin(target_hairpin_en)]} ({(min(target_hairpin_en)):.2f} kJ/mol). Algorithm may not converge."""
            OligamaWarning(mes, self)

        self.checked_candidates = set()  # Для предотвращения повторных проверок
        # Для поиска запускается Parallel. В этом режиме все потоки пишут в candidate_stats
        # Для синхронизации используется управляемый словарь 
        self.candidate_stats = {} 
        manager = Manager()
        self.candidate_stats = manager.dict()  # Создаем управляемый словарь
        self.failure_stats = manager.dict({  # Потокобезопасный словарь для статистики отказов
            "affinity_with_gene": 0,
            "gene_affinity_range": 0,
            "overall_complementary": 0,
            "consecutive_complementary": 0,
            "affinity_with_targets": 0,
            "hairpin_energy": 0
        })

        
        
    def find_olig_joblib(self):
        site_seq = self.gene[self.site_start:self.site_start + self.site_length]
        start_time = time.time()

        max_mutations = 10
        crossover_cycles = 2
        max_crossover_batch = 8

        iteration_count = 0

        found_oligos = []

        use_parallel = self.data.get("aff_predictor") == "Nupack"  # Условие для использования Parallel

        while time.time() - start_time < self.time_limit:

            iteration_count += 1

            candidates = [
                multiple_crossover(
                    [site_seq],
                    max_mutations=max_mutations,
                    crossover_cycles=crossover_cycles,
                    max_crossover_batch=max_crossover_batch
                ) for _ in range(self.batch_size)
            ]

            unique_candidates = [
                candidate for candidate in candidates if tuple(candidate) not in self.checked_candidates
            ]

            self.checked_candidates.update(tuple(candidate) for candidate in unique_candidates)


            if use_parallel:
                # Многопоточная обработка
                results = Parallel(n_jobs=4)(
                    delayed(self.check_candidate)(candidate) for candidate in unique_candidates
                )
            else:
                # Однопоточная обработка
                results = [self.check_candidate(candidate) for candidate in unique_candidates]

            for candidate, result in zip(unique_candidates, results):
                if result and candidate not in found_oligos:
                    found_oligos.append(candidate)
                    print(f"Valid candidate found: {candidate} at iteration {iteration_count}")

                    # Проверяем, что параметры кандидата сохраняются корректно
                    if candidate not in self.candidate_stats:
                        print(f"Error: Candidate {candidate} not found in candidate_stats.")
                        
                    if len(found_oligos) >= self.num_oligos:
                        #self.save_graphs(times, mutation_history, crossover_history, batch_history)
                        #self.save_stats()
                        return found_oligos



        print(f"Time limit exceeded: {self.time_limit} seconds. Iterations: {iteration_count}")
        OligamaWarning(f"Time limit exceeded: {self.time_limit} seconds. Iterations: {iteration_count}", self)
        #self.save_graphs(times, mutation_history, crossover_history, batch_history)
        #self.save_stats()
        return found_oligos


    def check_candidate(self, candidate):
        if not self.check_affinity_with_gene(candidate):
            self.failure_stats["affinity_with_gene"] += 1
            return False
        if not self.check_affinity_with_targets(candidate):
            self.failure_stats["affinity_with_targets"] += 1
            return False
        if not self.check_hairpin_energy(candidate):
            self.failure_stats["hairpin_energy"] += 1
            return False
        return True
      
      
    def check_affinity_with_gene(self, candidate):
        complementary_pairs = {"A": "T", "T": "A", "G": "C", "C": "G"}

        gene_site = self.gene[self.site_start:self.site_start + self.site_length]
        gene_aff = self.batch_affinity_predict([candidate], [gene_site])[0, 0]

        if not (float(self.gen_aff_low) <= gene_aff <= float(self.gen_aff_high)):
            self.failure_stats["gene_affinity_range"] += 1
            return False

        window_size = len(candidate)
        excluded_start = self.site_start - 1  # коррекция с человеческой интерпретации
        excluded_end = self.site_start + self.site_length

        # Создание реверсной последовательности
        reversed_candidate = "".join(reversed(candidate))

        best_overlap_count = 0
        best_max_consecutive = 0
        best_gene_segment = None

        for i in range(len(self.gene) - window_size + 1):
            if i + window_size > excluded_start and i < excluded_end:
                continue

            gene_segment = self.gene[i:i + window_size]
            overlap_count = 0
            max_consecutive = 0
            current_streak = 0

            for a, b in zip(reversed_candidate, gene_segment):
                if complementary_pairs.get(a) == b:
                    overlap_count += 1
                    current_streak += 1
                    max_consecutive = max(max_consecutive, current_streak)
                else:
                    current_streak = 0

            if overlap_count > best_overlap_count or (
                overlap_count == best_overlap_count and max_consecutive > best_max_consecutive
            ):
                best_overlap_count = overlap_count
                best_max_consecutive = max_consecutive
                best_gene_segment = gene_segment

            if overlap_count > self.max_overall_letters:
                self.failure_stats["overall_complementary"] += 1
                return False

            if max_consecutive > self.max_consecutive_letters:
                self.failure_stats["consecutive_complementary"] += 1
                return False

        # Сохранение параметров и участка гена
        self.candidate_stats[candidate] = {
            "overlap_count": best_overlap_count,
            "max_consecutive": best_max_consecutive,
            "gene_segment": best_gene_segment,
        }

        return True


    def save_stats(self):
        stats_path = self.output_folder / "gen_failure_stats.csv"

        failure_df = pd.DataFrame.from_dict(self.failure_stats, orient="index", columns=["failures"])
        failure_df.to_csv(stats_path)



    def check_affinity_with_targets(self, candidate):
        """
        Проверяет, соответствует ли кандидат условиям аффинности с target_seqs.
        """
        if not self.target_seqs.size:
            return True

        target_aff_arr = self.batch_affinity_predict(self.target_seqs, [candidate]).flatten()
        valid_low = (target_aff_arr > self.target_aff_low).all()
        valid_high = (target_aff_arr < self.target_aff_high).all()
        if valid_low and valid_high :
            print(f"affinity targets {candidate} OK")
        return valid_low and valid_high


    def check_hairpin_energy(self, candidate):
        """
        Проверяет, соответствует ли энергия шпильки кандидата пороговому значению.
        """
        hairpin_energy = self.hairpin_predictor.predict([candidate])[0]
        result =  hairpin_energy <= self.hairpin_en_thr
        print(f"hairpin energy {candidate}: {result}")
        return result

    def describe_oligos(self, seq_arr, toxls=True):
        """
        Создаёт описание олигонуклеотида, включая аффинности и энергии шпилек.
        """
        if not seq_arr or len(seq_arr) == 0:
            mes = "No oligos found. Please try less strict conditions."
            OligamaWarning(mes, self)
            return pd.DataFrame()

        # Приведение target_seqs к массиву
        if isinstance(seq_arr, str):
            seq_arr = [seq_arr]

        # Уникальные элементы в target_seqs
        seq_arr = list(set(seq_arr))

        try:
            # Создание списка всех последовательностей
            all_sequences = list([self.gene[self.site_start:self.site_start + self.site_length]]) + \
                            list(self.target_seqs) + \
                            list(seq_arr)

            # Расчёт полной матрицы аффинности
            affinity_matrix = self.batch_affinity_predict(all_sequences, all_sequences)

            # Формируем заголовки строк и столбцов
            labels = ["Gene Site"] + list(self.target_names) + list(seq_arr)

            # Преобразование матрицы аффинности в DataFrame
            affinity_df = pd.DataFrame(
                affinity_matrix,
                index=labels,
                columns=labels
            )

            # Добавление первой колонки для заголовков строк (дублирование индекса)
            affinity_df.insert(0, "Sequence", labels)



            # Сохранение результатов в Excel
            if toxls:

                # Создание таблицы параметров
                candidate_params = []
                for candidate in seq_arr:
                    params = self.candidate_stats.get(candidate, {
                        "overlap_count": None,
                        "max_consecutive": None,
                        "gene_segment": None,
                    })
                    candidate_params.append({
                        "Candidate": candidate,
                        "Overlap Count": params["overlap_count"],
                        "Max Consecutive": params["max_consecutive"],
                        "Gene Segment": params["gene_segment"],  # Добавляем участок гена
                    })

                params_df = pd.DataFrame(candidate_params)


                # Создание таблицы входных данных
                input_params = {
                    "Gene": self.gene,
                    "Site Start": self.site_start,
                    "Site Length": self.site_length,
                    "Max Consecutive Ratio (num)": f"{self.max_consecutive_ratio} ({self.max_consecutive_letters})",
                    "Max Overall Ratio (num)": f"{self.max_overall_ratio} ({self.max_overall_letters})",
                    "Gene Affinity Low": self.gen_aff_low,
                    "Gene Affinity High": self.gen_aff_high,
                    "Hairpin Energy Threshold": self.hairpin_en_thr,
                    "Metric": self.metric,
                    "Temperature (C)": self.celsius,
                    "Affinity Predictor": self.aff_predictor_name,  
                    "Hairpin Predictor": self.hairpin_predictor_name,  
                }
                input_params_df = pd.DataFrame(list(input_params.items()), columns=["Parameter", "Value"])

                # Объединение массивов в таблицу
                combined_data = {
                    "Target Names": self.target_names if len(self.target_names) > 0 else ["N/A"],
                    "Target Sequences (Orig)": self.target_seqs_orig if len(self.target_seqs_orig) > 0 else ["N/A"],
                    "Target Affinity Low": self.target_aff_low if len(self.target_aff_low) > 0 else ["N/A"],
                    "Target Affinity High": self.target_aff_high if len(self.target_aff_high) > 0 else ["N/A"],
                }

                # Преобразование объединенных данных в DataFrame
                max_length = max(len(combined_data[key]) for key in combined_data)  # Найти максимальную длину массива
                for key in combined_data:
                    # Заполняем массивы до максимальной длины пустыми строками для корректного отображения
                    combined_data[key] = list(combined_data[key]) + [""] * (max_length - len(combined_data[key]))

                target_df = pd.DataFrame(combined_data)


                output_path = self.output_folder / "Gene_Olig_Finder.xlsx"
                df_to_excel([affinity_df, params_df, input_params_df, target_df], ["Affinity Matrix", "Candidate Parameters", "Input Parameters", "Target Information"], output_path)
                
                #self.save_stats()
                
            return affinity_df

        except Exception as e:
            error_message = f"An error occurred while describing oligos: {e}"
            print(f"Error details: seq_arr={seq_arr}, target_seqs={self.target_seqs}")
            raise OligamaWarning(error_message, self)


    def batch_affinity_predict(self, seq1_arr, seq2_arr):
        master_seq_arr = self.arr_product(seq1_arr, seq2_arr)
        aff_arr = self.aff_predictor.predict(master_seq_arr[:, 0], master_seq_arr[:, 1], units=self.metric)
        aff_arr = aff_arr.reshape((len(seq1_arr), len(seq2_arr)))

        return aff_arr

    def arr_product(self, arr1, arr2):
        len1 = len(arr1)
        len2 = len(arr2)
        try:
            arr1 = np.array([arr1] * len2).reshape((len2, len1)).T
            arr2 = np.array([arr2] * len1).reshape((len1, len2))
        except ValueError as e:
            raise OligamaWarning(f"Error reshaping arrays: len1={len1}, len2={len2}, error={e}", self)

        return np.concatenate((arr1[:, :, None], arr2[:, :, None]), axis=2).reshape(-1, 2)
