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
    - Генерация кандидатов с мутациями антисенса до порогов target_aff_low/target_aff_high
    - Расчёт относительной аффинности через предикторы
    - Отсев по кросс-аффинностям и self-аффинности
    - Расчёт энергии и фильтр по Hairpin_energy_thr
    - Построение шаблона Pattern относительно комплемента к Target
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
        
        # Проверяем, что у нас только один target
        if len(self.target_seqs) != 1:
            OligamaWarning("AddOligsNew supports only one target sequence", self)
        if not self.target_seqs or not self.target_seqs[0]:
            OligamaWarning(f"Target sequence {self.target_names[0] if self.target_names else '[unknown]'} (index 0) is empty", self)
        self.target_seq = self.target_seqs[0] if self.target_seqs else ""
        self.target_name = self.target_names[0] if self.target_names else "target"

    def run(self, *, timeout_seconds: int = 300, min_required: int = 1) -> Dict[str, Any]:
        """
        Подбор кандидатов с ограничением по времени и минимальному количеству.

        Args:
            timeout_seconds (int): Максимум времени на подбор (по умолчанию 5 минут)
            min_required (int): Минимальное количество финальных кандидатов

        Returns:
            Dict[str, Any]: Словарь с результатами:
                - results: List[dict] - список финальных кандидатов
                - log: str - объединённый лог
        """
        # Основной цикл: пока не истечёт время или не наберём min_required
        deadline = time.time() + float(timeout_seconds)
        results: List[Dict[str, Any]] = []
        added_any = False

        while time.time() < deadline and len(results) < min_required:
            # Генерация кандидатов
            candidates = self._generate_candidates(
                self.target_seq,
                self.target_aff_low,
                self.target_aff_high,
                self.rounds
            )

            for candidate_seq in candidates:
                if time.time() >= deadline:
                    break

                # Лог заголовка теста
                self._log(f"\nTest {self.add_name}: {candidate_seq}")

                # Проверка кросс-аффинности с уже существующими в системе
                is_bad = False
                for other_seq in self._get_other_sequences():
                    rel_aff = self._compute_affinity(candidate_seq, other_seq)
                    self._log(f"Aff with {other_seq}: {self._format_float(rel_aff)}")
                    if rel_aff > self.bad_thr:
                        is_bad = True
                        break

                # Проверка self-аффинности
                if not is_bad:
                    self_aff = self._compute_affinity(candidate_seq, candidate_seq)
                    self._log(f"Aff with self - {candidate_seq}: {self._format_float(self_aff)}")
                    if self_aff > self.bad_thr:
                        is_bad = True

                if not is_bad:
                    # Энергия структуры кандидата
                    energy = self._compute_energy(candidate_seq)
                    if energy > self.Hairpin_energy_thr:
                        # Аффинность к таргету
                        rel_aff_to_target = self._compute_affinity(candidate_seq, self.target_seq)
                        
                        # Шаблон Pattern по совпадению с комплементарной таргету
                        pattern = self._build_pattern(candidate_seq, self._reverse_complement(self.target_seq))

                        results.append({
                            "name": self.add_name,
                            "seq": candidate_seq,
                            "Hairpin energy, kJ/mol": energy,
                            "affinity_to_target": rel_aff_to_target,
                            "pattern": pattern
                        })
                        added_any = True

            # Если кандидаты закончились и ничего не добавили — повторим генерацию при наличии времени
            if not added_any and time.time() < deadline:
                continue
            added_any = False

        return {
            "results": results,
            "log": "\n".join(self._log_lines)
        }

    def _generate_candidates(self, target_seq: str, min_thr: float, max_thr: float, rounds: int) -> List[str]:
        """Генерирует кандидатов на основе мутаций антисенса"""
        candidates = []
        asense_str = self._reverse_complement(target_seq)
        
        # Базовая аффинность целевого к антисенсу
        eq_rel = self._compute_affinity(target_seq, asense_str)
        
        letters_lower = ["a", "t", "g", "c"]
        letters_upper = ["A", "T", "G", "C"]
        length = len(asense_str)
        seen = set()

        for _ in range(rounds):
            str_mod = asense_str
            eq_rel_curr = eq_rel

            # Мутации, пока относительная аффинность > min_thr
            while eq_rel_curr > min_thr:
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

        return candidates

    def _compute_affinity(self, seq1: str, seq2: str) -> float:
        """Вычисляет относительную аффинность через предиктор"""
        return float(self.aff_predictor.predict([seq1], [seq2], units=self.metric)[0])

    def _compute_energy(self, seq: str) -> float:
        """Вычисляет энергию через предиктор"""
        return float(self.hairpin_predictor.predict([seq])[0])

    def _get_other_sequences(self) -> List[str]:
        """Возвращает список других последовательностей (кроме target)"""
        return [seq for seq, name in zip(self.target_seqs, self.target_names) 
                if name != self.target_name]

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

    def save_to_excel(self, df, filename="OligsFinderBash_results.xlsx"):
        df_to_excel(
            [df],
            ["Sheet1"],
            self.output_folder / filename,
            scientific_format_flag=False
        )        