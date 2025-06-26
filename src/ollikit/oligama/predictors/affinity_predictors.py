import pickle
import numpy as np
from keras import models
import nupack
from joblib import Parallel, delayed
import os
from oligama.utils import *
import logging


class CNN_Affinity_Predictor():

	def __init__(self,
				model_file = 'model_data/ensemble_CNN.keras',
				vect_layer_file = 'model_data/vect_layer.pkl', 
				celsius = 25):
		self.celsius = celsius
		if not os.path.exists(model_file):
			print("Файл модели не найден:", model_file)

		self.model = models.load_model(model_file, compile = False)
		self.unique_kmer_pairs = init_kmer_vocab(1)
		with open(vect_layer_file, 'rb') as vect_layer:
			self.vect_layer = pickle.load(vect_layer)

	def predictor_name(self) :
		return "Oligama"
    	
	def predict(self, seq1_arr, seq2_arr, units = 'fraction'):
		
		encoded = [encode_kmers(seq1, seq2, self.unique_kmer_pairs) for seq1, seq2 in zip(seq1_arr, seq2_arr)]
		encoded = self.vect_layer(encoded)
		fraction_arr = self.model(encoded)
		alpha_298_arr, alpha_310_arr = fraction_arr[:, 0], fraction_arr[:, 1]
		return self.convert_results(alpha_298_arr, alpha_310_arr, units, self.celsius)
		  
	def convert_results(self, alpha_298_arr, alpha_310_arr, units, celsius):
		results = np.empty_like(alpha_298_arr)
		if units == 'fraction' and celsius == 25:
			results = alpha_298_arr
		elif units == 'fraction' and celsius == 37:
			results = alpha_310_arr
		else:

			k_298_array = 1 / alpha_298_arr + alpha_298_arr - 2
			k_310_array = 1 / alpha_310_arr + alpha_310_arr - 2


			G_298_array = (np.log(k_298_array) - 6 * np.log(10)) * 8.31 * 298 * 10**(-3)
			G_310_array = (np.log(k_310_array) - 6 * np.log(10)) * 8.31 * 298 * 10**(-3)


			S_array = (G_310_array - G_298_array) / (298 - 310)
			H_array = (298 * G_310_array - 310 * G_298_array) / (298 - 310)
			T = celsius + 273
			G_array = H_array - T * S_array


			k_array = np.exp(G_array * 10**3 / 8.31 / 298) * 10**6
			alpha_array = 0.5 * ((2 + k_array) - np.sqrt(k_array**2 + 4 * k_array))


			if units == 'Kd' and celsius == 25:
				results = k_298_array/1e6
			elif units == 'Kd' and celsius == 37:
				results = k_310_array/1e6
			elif units == 'gibbs' and celsius == 25:
				results = G_298_array*1e3
			elif units == 'gibbs' and celsius == 37:
				results = G_310_array*1e3
			elif units == 'gibbs':
				results = G_array*1e3
			elif units == 'Kd':
				results = k_array/1e6
			elif units == 'fraction':
				results = alpha_array

		return np.array(results)


class Nupack_Affinity_Predictor():
	
	def __init__(self, celsius = 25, material = 'dna'):
		import nupack
		self.celsius = celsius
		self.model = nupack.Model(material=material, celsius = celsius)
		logging.info(f"Инициализирован Nupack_Affinity_Predictor с параметрами: celsius={celsius}, material={material}")

	def predictor_name(self) :
		return "Nupack"

	def simple_aff(self, seq1, seq2,  units = 'fraction'):
		logging.info(f"Расчет аффинности для пары последовательностей:\nseq1: {seq1}\nseq2: {seq2}")
		s1 = nupack.Strand(seq1, name='s1')
		s2 = nupack.Strand(seq2, name='s2')

		s1_s2 = nupack.Complex([s1, s2])
		complexes = nupack.SetSpec(max_size=2)

		complex_set = nupack.ComplexSet(strands=[s1,s2], complexes=complexes)
		logging.info("Начало анализа комплексов...")
		result = nupack.complex_analysis(complexes=complex_set, model=self.model, compute=['pfunc', 'mfe'])
		logging.info("Анализ комплексов завершен")
		
		Kw = nupack.constants.water_molarity(self.model.temperature)
		logging.info(f"Константа воды (Kw): {Kw}")

		pfunc_s1_s2 = result[s1_s2].pfunc
		if pfunc_s1_s2 == 0:
			logging.error(f"Partition function для {s1_s2} равна нулю. Невозможно вычислить аффинность.")
			raise ValueError(f"Partition function for {s1_s2} is zero. Cannot compute affinity.")
		
		pfunc_s1 = result['(s1)'].pfunc
		pfunc_s2 = result['(s2)'].pfunc
		logging.info(f"Partition functions:\npfunc_s1: {pfunc_s1}\npfunc_s2: {pfunc_s2}\npfunc_s1_s2: {pfunc_s1_s2}")

		k = Kw * float(pfunc_s1 * pfunc_s2 / pfunc_s1_s2)
		logging.info(f"Константа равновесия (k): {k}")

		if units == 'fraction':
			k *= 10**6
			result = 0.5*((2+k) - np.sqrt(k**2 + 4*k))
			logging.info(f"Результат в единицах fraction: {result}")
			return result
		elif units == 'Kd':
			logging.info(f"Результат в единицах Kd: {k}")
			return k
		elif units == 'gibbs':
			gibbs = 8.31*(273+self.celsius)*np.log(k)
			logging.info(f"Результат в единицах gibbs: {gibbs}")
			return gibbs


	def predict(self, seq1_arr, seq2_arr, units='fraction', n_jobs=-1):
		logging.info(f"Начало предсказания аффинности для {len(seq1_arr)} пар последовательностей")
		logging.info(f"Используется {n_jobs} параллельных процессов")
		
		def safe_simple_aff(seq1, seq2):
			try:
				result = self.simple_aff(seq1, seq2, units)
				logging.info(f"Успешно рассчитана аффинность для пары:\n{seq1}\n{seq2}\nРезультат: {result}")
				return result
			except ValueError as e:
				logging.error(f"Ошибка при расчете аффинности для пары ({seq1}, {seq2}): {e}")
				return np.nan

		results = np.array(
			Parallel(n_jobs=n_jobs)(
				delayed(safe_simple_aff)(seq1, seq2) for seq1, seq2 in zip(seq1_arr, seq2_arr)
			)
		)
		
		nan_count = np.isnan(results).sum()
		if nan_count > 0:
			logging.warning(f"Обнаружено {nan_count} NaN значений в результатах")
		
		logging.info("Предсказание аффинности завершено")
		return results


