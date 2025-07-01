import hashlib
import json
from ollikit import Concentration_Solver  # Импортируем только нужный класс
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.optimize import minimize
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import os

count_seq = 8  # Количество олигонуклеотидов

# Входные данные
input_template = {
    "target_names": [],
    "sequence_types": ["DNA"] * count_seq,  # Максимальное количество олигонуклеотидов, можно менять
    "total_conc": ["1e-6"] * count_seq,  # Максимальное количество олигонуклеотидов, можно менять
    "input_seq_name": "s1",
    "output_seq_name": "s2",
    "input_seq_bounds": ["1e-9", "1e-6"],
    "n_points": 20,
    "aff_predictor": "Oligama",  # Меняем на Oligama, чтобы не требовать Nupack
    "hairpin_predictor": "Seqfold",
    "metric": "fraction",
    "celsius": 25
}


# Функция для генерации MD5 хеша из последовательностей
def generate_md5_from_sequences(sequences):
    """
    Генерирует MD5 хеш из списка последовательностей.
    """
    hash_object = hashlib.md5(''.join(sequences).encode())
    return hash_object.hexdigest()


# Функция для сохранения input_data в файл
def save_input_data(input_data, output_folder, hash_value):
    """
    Сохраняет input_data в текстовый файл с названием, основанным на MD5 хеше.
    """
    file_path = output_folder / Path(f"{hash_value}.json")
    with open(file_path, 'w') as f:
        json.dump(input_data, f, indent=4)


# Функция для генерации случайной DNA последовательности
def generate_random_dna_sequence(length=20):
    """
    Генерирует случайную DNA последовательность длины `length`.
    Состав: случайный выбор из A, T, G, C.
    """
    return ''.join(np.random.choice(['A', 'T', 'G', 'C'], size=length))


# Функция для запуска концентрационного решения
def run_concentration_solver(input_file, output_folder):
    # Запуск предсказания
    solver = Concentration_Solver(input_file, output_folder)
    result, inp_df = solver.predict_only()  # Результат возвращается как два DataFrame
    return result, inp_df


# Функция для построения графика
def plot_concentration(result, output_folder, input_data):
    # Извлекаем концентрации из результата
    inp_conc = result[f"{result.columns[0]}"]  # Входная концентрация
    out_conc = result[f"{result.columns[1]}"]  # Выходная концентрация

    # Оценка, насколько зависимость похожа на синусоиду
    # sinusoidal_error = evaluate_sinusoidal_fit(out_conc)
    sinusoidal_error, optimal_phase = evaluate_sinusoidal_fit_with_phase(out_conc)
    # Если ошибка синусоидальной подгонки слишком велика, пропускаем сохранение
    if sinusoidal_error > 18:
        #print(f"Ошибка синусоидальной подгонки слишком велика: {sinusoidal_error}, пропускаем сохранение.")
        return
    
    # Строим график
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(inp_conc, out_conc, marker='o', linestyle='-', color='b', label="Концентрация выхода")
    ax.scatter(inp_conc, out_conc)
    
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(f"{result.columns[0]} concentration, M")
    ax.set_ylabel(f"{result.columns[1]} concentration, M")
    ax.legend()
    
    # Генерация уникального хеша для имени файла
    hash_value = generate_md5_from_sequences(input_data['target_seqs'])

    # Создание имени файла с ошибкой синусоидальности
    file_name = f"{sinusoidal_error:.4f}_{hash_value}"

    # Сохраняем график в папку output с уникальным именем
    plt.savefig(output_folder / Path(f"{file_name}.png"))
    plt.close()

    # Сохраняем input_data в файл с тем же именем, что и у графика
    save_input_data(input_data, output_folder, file_name)


# Функция для оценки, насколько зависимость похожа на синусоиду
def evaluate_sinusoidal_fit(concentrations):
    """
    Оценка, насколько зависимость концентрации приближена к синусоиде.
    Мы будем использовать минимальную ошибку при подгонке синусоиды.
    """
    # Создаём синусоиду с той же длиной
    x = np.arange(len(concentrations))
    y_sin = np.sin(x)  # Периодическая функция (синусоида)

    # Нормализуем данные (для сравнения)
    concentrations = np.array(concentrations) / np.max(concentrations)
    y_sin = y_sin / np.max(np.abs(y_sin))

    # Считаем ошибку подгона
    error = np.sum((concentrations - y_sin) ** 2)
    return error



# Функция для оценки, насколько зависимость похожа на синусоиду с учётом фазы
def evaluate_sinusoidal_fit_with_phase(concentrations):
    """
    Оценка, насколько зависимость концентрации приближена к синусоиде.
    Мы будем использовать минимизацию ошибки с учётом фазы синусоиды.
    """
    # Функция ошибки для минимизации
    def error_function(phase, concentrations):
        x = np.arange(len(concentrations))
        # Создаем синусоиду с учётом фазы
        y_sin = np.sin(x + phase)
        # Нормализуем данные (для сравнения)
        concentrations = np.array(concentrations) / np.max(concentrations)
        y_sin = y_sin / np.max(np.abs(y_sin))
        # Возвращаем сумму квадратов ошибок (чем меньше, тем лучше)
        return np.sum((concentrations - y_sin) ** 2)

    # Начальное значение фазы
    initial_phase = 0.0
    # Минимизация ошибки с учётом фазы
    result = minimize(error_function, initial_phase, args=(concentrations,), bounds=[(-np.pi, np.pi)])

    # Результат минимизации - это оптимальная фаза
    optimal_phase = result.x[0]

    # Считаем ошибку подгона с найденной фазой
    error = result.fun
    return error, optimal_phase

# Пример использования
# concentrations = np.random.rand(100)  # Пример случайных данных
# error, optimal_phase = evaluate_sinusoidal_fit_with_phase(concentrations)

# print(f"Ошибка синусоидальной подгонки: {error}")
# print(f"Оптимальная фаза: {optimal_phase}")


# Функция для перебора олигонуклеотидов в параллельных потоках
def generate_oligo_combinations(input_template, output_folder, max_oligos=8, num_threads=10):
    # Параллельное выполнение задач с использованием ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        futures = []
        
        # Запускаем генерацию олигонуклеотидов в многозадачном режиме
        for _ in range(max_oligos):
            futures.append(executor.submit(process_oligo_set, input_template, output_folder))
        
        # Ожидаем завершения всех задач
        for future in futures:
            future.result()


# Функция для обработки одного набора олигонуклеотидов
def process_oligo_set(input_template, output_folder):
    input_data = input_template.copy()
    while True:


        # Генерация новых последовательностей и имен
        input_data['target_seqs'] = [generate_random_dna_sequence(20) for _ in range(count_seq)]  # np.random.randint(15, 25)
        input_data['target_names'] = [f"s{i+1}" for i in range(count_seq)]
        
        # Запуск решения для текущего набора олигонуклеотидов
        result, inp_df = run_concentration_solver(input_data, output_folder)

        # Построение и сохранение графика
        plot_concentration(result, output_folder, input_data)


# Пример использования
output_folder = Path("test_output/conc_solver_visual")  # Папка для сохранения графиков
os.makedirs(output_folder, exist_ok=True)  # Создаём папку, если её нет

generate_oligo_combinations(input_template, output_folder, max_oligos=count_seq, num_threads=11)
