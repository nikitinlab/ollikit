# ollikit

[![PyPI - Version](https://img.shields.io/pypi/v/ollikit.svg)](https://pypi.org/project/ollikit)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/ollikit.svg)](https://pypi.org/project/ollikit)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE.txt)

---

**ollikit** — современный Python-пакет для биоинформатики, анализа олигонуклеотидов, расчёта аффинности, работы с физико-химическими и нейросетевыми предикторами, а также решения сложных задач концентраций и поиска оптимальных последовательностей.

---

## 📦 Возможности

- Расчёт матриц аффинности между олигонуклеотидами (Affinity_Matrix)
- Предсказание энергии шпилек (hairpin) и комплексообразования
- Решение систем концентраций (в том числе с несколькими входами/выходами)
- Поиск оптимальных олигонуклеотидов по заданным критериям (Olig_Finder, Gene_Olig_Finder)
- Поддержка различных предикторов (Nupack, Seqfold, Oligama и др.)
- Генерация и анализ случайных последовательностей (Random_Initializer)
- Работа с microRNA и генными олигонуклеотидами (MicroRNA_Finder)
- Удобное сохранение результатов (Excel, графики)
- Современная архитектура, покрытие тестами, расширяемость

---

## 🚀 Быстрый старт

```python
from ollikit import Affinity_Matrix, Complex_Solver, Olig_Finder

input_data = {
    "target_seqs": ["ACTGCTAGAGATTTTCCACAT", "TGTGGAAAATCTCTAGCAGTT"],
    "target_names": ["s1", "s2"],
    "metric": "fraction",
    "celsius": 25,
    "aff_predictor": "Nupack",
    "hairpin_predictor": "Seqfold"
}
output_folder = "results"
aff = Affinity_Matrix(input_data, output_folder)
aff_matr = aff.predict()
aff.heatmaps_to_png(aff_matr)
aff.matrix_to_excel(aff_matr)
```

---

## 🧬 Основные классы и функции

- `Affinity_Matrix` — построение и визуализация матриц аффинности
- `Complex_Solver` — расчёт равновесных концентраций в сложных смесях
- `Olig_Finder` — поиск оптимальных олигонуклеотидов
- `Primary_Concentration_Solver`, `Concentration_Solver` — анализ зависимости выхода от входа
- `Gene_Olig_Finder`, `MicroRNA_Finder` — работа с генными и микроРНК-мишенями
- `Random_Initializer` — генерация случайных последовательностей
- `Complement` — операции с комплементарностью
- `Nupack_Affinity_Predictor`, `Nupack_Hairpin_Predictor` — физико-химические предикторы
- `triangle_to_matr`, `detect_unit`, `convert`, `fpi_step`, `find_equilibrium_conc` — утилиты и вспомогательные функции

---


## 🛠️ Установка

### С PyPI

```bash
pip install https://github.com/nikitinlab/ollikit/raw/refs/heads/master/dist/ollikit-0.1.1.tar.gz
# pip install ollikit

```

### Из исходников

```bash
git clone https://github.com/nikitinlab/ollikit.git
cd ollikit
pip install -e .
```

---

## ⚡ Зависимости

- numpy, pandas, matplotlib, seaborn, scipy, biopython, numba, joblib, xlsxwriter, nupack, seqfold, tensorflow, keras и др.

Для отключения warnings  используйте fork seqfold:
pip install https://github.com/nikitinlab/ollikit/raw/refs/heads/master/dist/seqfold-0.7.18.post1-py3-none-any.whl 
---

## 🧪 Тестирование

```bash
pytest
```
или для конкретной папки:
```bash
pytest ollikit/tests/
```


## 💡 Вклад и поддержка

- Pull requests и баг-репорты приветствуются!
- Для связи: [email](mailto:725156@gmail.com)
- [Issues](https://github.com/nikitinlab/ollikit/issues)

---

## 📄 Лицензия

MIT License

---

**ollikit** — универсальный инструмент для анализа олигонуклеотидов и биоинформатики!
