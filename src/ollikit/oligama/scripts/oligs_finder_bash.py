import numpy as np
import pandas as pd
import time
from pathlib import Path
from typing import List, Dict, Any

from ..data_loaders import AddOligsNew_Dataloader
from ..utils import df_to_excel
from ..exceptions import OligamaWarning, OligamaException

class AddOligsNew(AddOligsNew_Dataloader):

    """
    Класс для подбора новых олигонуклеотидов на основе порогов аффинности и энергии.
    
    Основные возможности:
    - Генерация кандидатов с мутациями антисенса до порогов affinity_low/affinity_high
    - Расчёт относительной аффинности через предикторы
    - Отсев по кросс-аффинностям и self-аффинности
    - Расчёт энергии и фильтр по Hairpin_energy_thr
    - Построение шаблона Pattern относительно комплемента к Target
    - Target определяется по полю target_name из списка target_names
    """

    def __init__(self, input_data, output_folder):
        """
        Инициализирует объект AddOligsNew.

        Args:
            input_data (str или dict): Путь к файлу JSON или словарь с входными данными
            output_folder (str): Путь к папке для выходных данных
        """
        super().__init__(input_data, output_folder)
        
        # Буферы в памяти
        self._log_lines: List[str] = []
        self._rand_candidates: List[str] = []
        
        # Найти target по target_name среди всех последовательностей
        if not hasattr(self, 'target_seqs') or len(self.target_seqs) == 0:
            OligamaWarning("No sequences provided", self)
        
        # Преобразуем в списки если это numpy arrays
        target_seqs_list = list(self.target_seqs) if hasattr(self.target_seqs, '__iter__') else [self.target_seqs]
        target_names_list = list(self.target_names) if hasattr(self.target_names, '__iter__') else [self.target_names]
        
        # Ищем индекс target'а по target_name
        target_index = None
        try:
            target_index = target_names_list.index(self.target_name)
        except ValueError:
            OligamaWarning(f"Target name '{self.target_name}' not found in target_names: {self.target_names}", self)
            # Используем первую последовательность как fallback
            target_index = 0
        
        if target_index is None or target_index >= len(target_seqs_list):
            OligamaWarning(f"Invalid target index {target_index}", self)
            target_index = 0
        
        # Устанавливаем target последовательность
        self.target_seq = target_seqs_list[target_index] if target_seqs_list else ""
        
        # Контрольные последовательности (все кроме target'а)
        self.control_seqs = [seq for i, seq in enumerate(target_seqs_list) if i != target_index]
        self.control_names = [name for i, name in enumerate(target_names_list) if i != target_index]

    def run(self, *, timeout_seconds: int = None, min_required: int = None) -> Dict[str, Any]:
        """
        Подбор кандидатов с ограничением по времени и минимальному количеству.

        Args:
            timeout_seconds (int): Максимум времени на подбор (секунды, если None — берётся self.timeout из формы)
            min_required (int): Минимальное количество финальных кандидатов (если None — берётся self.num_oligos)

        Returns:
            Dict[str, Any]: Словарь с результатами:
                - results: List[dict] - список финальных кандидатов
                - log: str - объединённый лог
        """
        import datetime
        log_path = self.output_folder / "Log.txt"
        def write_log(msg):
            with open(log_path, "a") as f:
                f.write(f"[{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {msg}\n")

        # Логируем стартовые параметры
        if min_required is None:
            min_required = self.num_oligos
        write_log(f"START AddOligsNew.run: timeout={timeout_seconds if timeout_seconds is not None else self.timeout}, min_required={min_required} (from num_oligos={self.num_oligos}), target_name={self.target_name}, target_seq={self.target_seq}, control_seqs_count={len(self.control_seqs)}, affinity_low={self.affinity_low}, affinity_high={self.affinity_high}, bad_thr={self.bad_thr}, rounds={self.rounds}, Hairpin_energy_thr={self.Hairpin_energy_thr}")

        if timeout_seconds is None:
            timeout_seconds = float(self.timeout)
        deadline = time.time() + float(timeout_seconds)
        # Устанавливаем deadline как атрибут экземпляра для доступа из _generate_candidates
        self._deadline = deadline
        results: List[Dict[str, Any]] = []
        added_any = False
        iteration = 0

        while time.time() < deadline and len(results) < min_required:
            iteration += 1
            write_log(f"ITERATION {iteration}: candidates generation started")
            candidates = self._generate_candidates(
                self.target_seq,
                self.affinity_low,
                self.affinity_high,
                self.rounds
            )
            write_log(f"ITERATION {iteration}: candidates generated: {len(candidates)}")

            found_this_iter = 0
            for candidate_seq in candidates:
                if time.time() >= deadline:
                    write_log(f"ITERATION {iteration}: deadline reached, breaking candidate loop")
                    break

                # Проверка кросс-аффинности с контрольными последовательностями
                is_bad = False
                bad_reason = ""
                for i, other_seq in enumerate(self._get_other_sequences()):
                    rel_aff = self._compute_affinity(candidate_seq, other_seq)
                    if rel_aff > self.bad_thr:
                        is_bad = True
                        control_name = self.control_names[i] if i < len(self.control_names) else f"control_{i}"
                        bad_reason = f"cross-affinity {rel_aff:.6f} > {self.bad_thr} with control {control_name}: {other_seq}"
                        break

                # Проверка self-аффинности
                if not is_bad:
                    self_aff = self._compute_affinity(candidate_seq, candidate_seq)
                    if self_aff > self.bad_thr:
                        is_bad = True
                        bad_reason = f"self-affinity {self_aff:.6f} > {self.bad_thr}"

                if not is_bad:
                    # Энергия структуры кандидата
                    try:
                        energy = self._compute_energy(candidate_seq)
                        if energy is None:
                            energy = 0.0
                    except Exception:
                        energy = 0.0                    
                    if energy > self.Hairpin_energy_thr:
                        rel_aff_to_target = self._compute_affinity(candidate_seq, self.target_seq)
                        pattern = self._build_pattern(candidate_seq, self._reverse_complement(self.target_seq))
                        results.append({
                            "Name": self.add_name,
                            "Sequence": candidate_seq,
                            "Hairpin, kJ/mol": energy,
                            "Affinity": rel_aff_to_target,
                            "Pattern": pattern
                        })
                        found_this_iter += 1
                        write_log(f"ITERATION {iteration}: FOUND candidate: {candidate_seq}, Hairpin={energy}, Affinity={rel_aff_to_target}, Pattern={pattern}")
                        added_any = True
                else:
                    write_log(f"ITERATION {iteration}: REJECTED candidate: {candidate_seq}, reason: {bad_reason}")

            write_log(f"ITERATION {iteration}: found {found_this_iter} candidates, total found: {len(results)}")

            # Если кандидаты закончились и ничего не добавили — повторим генерацию при наличии времени
            if not added_any and time.time() < deadline:
                write_log(f"ITERATION {iteration}: no candidates found, retrying generation")
                continue
            added_any = False

        write_log(f"FINISH AddOligsNew.run: total found={len(results)}, time elapsed={int(time.time() + float(timeout_seconds) - deadline)}s")
        return {
            "results": results,
            "log": "\n".join(self._log_lines)
        }

    def _generate_candidates(self, target_seq: str, min_thr: float, max_thr: float, rounds: int) -> List[str]:
        """Генерирует кандидатов на основе мутаций антисенса"""
        import datetime
        log_path = self.output_folder / "Log.txt"
        def write_log(msg):
            with open(log_path, "a") as f:
                f.write(f"[{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {msg}\n")
        
        candidates = []
        asense_str = self._reverse_complement(target_seq)
        
        # Базовая аффинность целевого к антисенсу
        eq_rel = self._compute_affinity(target_seq, asense_str)
        write_log(f"_generate_candidates: asense={asense_str}, initial_eq_rel={eq_rel}, min_thr={min_thr}, max_thr={max_thr}")
        
        letters_lower = ["a", "t", "g", "c"]
        letters_upper = ["A", "T", "G", "C"]
        length = len(asense_str)
        seen = set()

        for round_num in range(rounds):
            str_mod = asense_str
            eq_rel_curr = eq_rel
            mutations_count = 0

            # Мутации, пока относительная аффинность > min_thr
            while eq_rel_curr > min_thr:
                mutations_count += 1
                
                # Логируем каждые 20 циклов
                if mutations_count % 20 == 0:
                    write_log(f"Round {round_num+1}/{rounds}, mutation {mutations_count}: eq_rel_curr={eq_rel_curr}, str_mod={str_mod}")

                # Проверяем timeout каждые 100 мутаций
                if mutations_count % 100 == 0:
                    deadline = getattr(self, '_deadline', None)
                    if deadline is not None and time.time() >= deadline:
                        write_log(f"Round {round_num+1}: timeout reached at mutation {mutations_count}, breaking mutation loop")
                        break


                idx = np.random.randint(length)
                base_idx = np.random.randint(4)
                pick_upper = letters_upper[base_idx]
                pick_lower = letters_lower[base_idx]

                # Если одинаково с исходной буквой -> верхний регистр, иначе нижний
                repl = pick_upper if pick_upper == str_mod[idx:idx+1] else pick_lower
                str_mod = str_mod[:idx] + repl + str_mod[idx+1:]

                # Считаем аффинность для текущей версии
                eq_rel_curr = self._compute_affinity(target_seq, str_mod)

                # Если попали в окно (min_thr, max_thr) -> сохраняем
                if min_thr < eq_rel_curr < max_thr:
                    str_cap = str_mod.upper()
                    if str_cap not in seen:
                        candidates.append(str_cap)
                        seen.add(str_cap)
                        write_log(f"Round {round_num+1}: found candidate {str_cap} with affinity {eq_rel_curr}")
                
                # Защита от бесконечного цикла
                if mutations_count > 1000:
                    write_log(f"Round {round_num+1}: breaking after 1000 mutations, current affinity={eq_rel_curr}")
                    break
            
            write_log(f"Round {round_num+1} completed: {mutations_count} mutations, candidates_so_far={len(candidates)}")

        write_log(f"_generate_candidates completed: total candidates={len(candidates)}")
        return candidates

    def _compute_affinity(self, seq1: str, seq2: str) -> float:
        """Вычисляет относительную аффинность через предиктор"""
        result = float(self.aff_predictor.predict([seq1], [seq2], units=self.metric)[0])
        # Диагностика: логируем первые несколько вызовов
        if not hasattr(self, '_affinity_calls'):
            self._affinity_calls = 0
        self._affinity_calls += 1
        if self._affinity_calls <= 5:
            import datetime
            log_path = self.output_folder / "Log.txt"
            with open(log_path, "a") as f:
                f.write(f"[{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] _compute_affinity: seq1={seq1}, seq2={seq2}, result={result}\n")
        return result

    def _compute_energy(self, seq: str) -> float:
        """Вычисляет энергию через предиктор"""
        return float(self.hairpin_predictor.predict([seq])[0])

    def _get_other_sequences(self) -> List[str]:
        """Возвращает список контрольных последовательностей (все кроме target)"""
        return self.control_seqs

    @staticmethod
    def _reverse_complement(seq: str) -> str:
        """Возвращает обратную комплементарную последовательность"""
        comp = {"A": "T", "T": "A", "G": "C", "C": "G", 
                "a": "t", "t": "a", "g": "c", "c": "g"}
        return "".join(comp.get(base, base) for base in reversed(seq))

    @staticmethod
    def _build_pattern(word1: str, word2: str, wrong_char: str = "x") -> str:
        """Строит шаблон совпадения последовательностей"""
        pattern_chars = []
        len2 = len(word2)
        for j in range(len2):
            c1 = word1[j] if j < len(word1) else ""
            c2 = word2[j]
            pattern_chars.append(c1 if c1 == c2 and c1 != "" else wrong_char)
        return "".join(pattern_chars)

    @staticmethod
    def _format_float(value: float) -> str:
        """Форматирует число с плавающей точкой"""
        return f"{value:g}"

    def _log(self, text: str) -> None:
        """Добавляет текст в лог"""
        self._log_lines.append(text)

    def save_to_excel(self, df, filename="OligsFinder2_results.xlsx"):
        df_to_excel(
            [df],
            ["Sheet1"],
            self.output_folder / filename,
            scientific_format_flag=False
        )        