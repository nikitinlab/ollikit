import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os


from ..scripts.affinity_matrix import Affinity_Matrix
from ..utils import df_to_excel
from ..data_loaders import *
from ..complex_concentration_minimal import find_eq_conc
from ..predictors.hairpin_predictors import *

class Complex_Solver(Complex_Solver_Dataloader):
	def __init__(self, input_file, output_folder):
		super().__init__(input_file, output_folder)

		# Создание матрицы аффинности
		Aff_matr = Affinity_Matrix(input_file, output_folder)
		self.aff_matr = Aff_matr.predict(visualize=False, units='gibbs')

		n = len(self.target_seqs)
		self.linear_g = np.concatenate((np.zeros(n), self.aff_matr[np.triu_indices(n, 1)]))
  
		if self.hairpin_predictor_name == 'Seqfold':
			self.hairpin_predictor = Seqfold_Hairpin_Predictor(celsius = self.celsius)

	def predict(self, toxls=True):
		# Получаем концентрации
		complex_concs = find_eq_conc(self.linear_g, self.total_conc, self.celsius)

		# Число последовательностей
		n = len(self.target_seqs)

		# Разделение на индивидуальные концентрации и комплексы
		individual_concs = complex_concs[:n]
		complex_concs_only = complex_concs[n:]

		# Надо посчитать энергию шпилек и выдать предупрежение
		hairpin_cond = np.array(self.hairpin_predictor.predict(self.target_seqs))  
  
		hairpin_normal_limit = -0.9 # ЛИМИТ НОРМАЛЬНОСТИ ДЛЯ HAIRPIN (символы !!! меняют цвет столбца на красный)
		hairpin_with_flags = [
			f"{cond} !!!" if cond < hairpin_normal_limit else cond for cond in hairpin_cond
        ]
		individual_data = {
            'Type': 'Individual',
            'Name': self.target_names,  # Используем заданные имена последовательностей
            'Concentration': individual_concs,
            'Hairpin': hairpin_with_flags
        }
		complexes_data = {
			'Type': 'Complex',
			'Name': [f"{self.target_names[i]}-{self.target_names[j]}" 
					for i, j in zip(*np.triu_indices(n, 1))],
			'Concentration': complex_concs_only,
			'Hairpin':  [""] * len(complex_concs_only) 
		}

		data = pd.concat([pd.DataFrame(individual_data), pd.DataFrame(complexes_data)])
		
  
		# Экспорт таблицы в Excel
		if toxls:
			excel_path = os.path.join(self.output_folder, "data.xlsx")
			df_to_excel([data], ['Complex'], excel_path)


		# Импорт библиотек для стиля
		from matplotlib.colors import LinearSegmentedColormap

		# Горизонтальный столбчатый график
		plt.figure(figsize=(12, 8))

		# Данные для графика
		data = data.iloc[::-1]
		names = data['Name']
		concentrations = data['Concentration']
		colors = ['red' if (row.Type == 'Individual' and row.Hairpin != '' and 
                            isinstance(row.Hairpin, str) and "!!!" in row.Hairpin) 
                  else 'darkblue' if row.Type == 'Individual' else 'maroon'
                  for _, row in data.iterrows()]

		# Построение горизонтального столбчатого графика
		bars = plt.barh(names, concentrations, color=colors, alpha=0.9, edgecolor='black')

		# Градиентные заливки
		for i, bar in enumerate(bars):
			bar_color = colors[i]
			cmap = LinearSegmentedColormap.from_list("gradient", [bar_color, "white"])
			bar.set_facecolor(cmap(0.5))  # Применяем градиент
			bar.set_alpha(0.8)  # Делаем столбики слегка прозрачными

		# Полосы для лучшей читаемости
		for i in range(0, len(names), 2):
			plt.axhspan(i - 0.5, i + 0.5, facecolor='lightgrey', alpha=0.2)

		# Половина размаха графика
		range_half = (max(concentrations) - min(concentrations)) / 2

		# Подписи значений на каждом столбике
		for i, bar in enumerate(bars):
			width = bar.get_width()  # Длина столбика
			name = data.iloc[i]['Name']
			hairpin_value = data.loc[data['Name'] == name, 'Hairpin'].values
			# Проверяем, есть ли значение и соответствует ли условию
			if hairpin_value.size > 0 and colors[i] == 'red':
				value = hairpin_value[0]
				if isinstance(value, str) and "!!!" in value:
					try:
						value = float(value.split()[0])  # Преобразуем число перед "!!!" в float
					except ValueError:
						value = None  # На случай, если значение не преобразуется
				Gh = f"Gh={value:.2f}" if isinstance(value, (int, float)) else ''
			else:
				Gh = ''
    
			if width > range_half:
				# Если значение большое, пишем внутри столбика
				plt.text(width - max(concentrations) * 0.01,  # Чуть левее конца столбика
						bar.get_y() + bar.get_height() / 2,  # Центр столбика
						f"{width:.2e}  {Gh}",  # Значение концентрации
						va='center', ha='right', fontsize=14, weight='normal', color='black')
			else:
				
				# Если значение небольшое, пишем снаружи
				plt.text(width + max(concentrations) * 0.01,  # Чуть правее столбика
						bar.get_y() + bar.get_height() / 2,  # Центр столбика
						f"{width:.2e} {Gh}",  # Значение концентрации
						va='center', ha='left', fontsize=14, weight='normal', color='black')


		# Настройка осей
		plt.xlabel('Concentration', fontsize=14, weight='bold', color='black')
		plt.ylabel('Sequences and Complexes', fontsize=14, weight='bold', color='black')
		plt.title('Complex Solver Results', fontsize=16, weight='bold', color='black')

		# Увеличение расстояния между метками на оси Y
		plt.yticks(fontsize=12, rotation=0)

		# Сетка
		plt.grid(axis='x', linestyle='--', alpha=0.4, color='grey')

		# Сохранение графика
		plot_path = os.path.join(self.output_folder, "bars.png")
		plt.tight_layout()
		plt.savefig(plot_path)
		plt.close()

		# Возвращаем таблицу для дальнейшего использования
		return data

