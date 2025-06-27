
import pandas as pd
from pathlib import Path
import json

from ..utils import df_to_excel, random_seq, GC_content
from ..exceptions import OligamaException, OligamaWarning
from ..predictors.hairpin_predictors import *


class Random_Initializer():
	def __init__(self, input_data, output_folder):
		self.output_folder = Path(output_folder)
		self.log_file = output_folder/Path("Log.txt")
  
		if isinstance(input_data, str):
            # Загружаем данные из файла
			with open(input_data) as f:
				data = json.load(f)
		elif isinstance(input_data, dict):
            # Используем переданные JSON-данные
			data = input_data
		else:
			raise OligamaException("input_data должен быть либо путем к файлу, либо JSON-объектом", self)
   
		self.GC_content_bounds = [float(x) for x in data["GC_content_bounds"]]
		self.num_oligos = data["num_oligos"]
		self.oligo_length = data["oligo_length"]
		self.celsius = data["celsius"]
		if self.celsius > 100 or self.celsius < 0:
			raise OligamaException("Temperature must be < 100 and > 0", self)
		self.hairpin_predictor_name = data["hairpin_predictor"]
		self.hairpin_thr = float(data["hairpin_thr"])
		self.output_format = data["output_format"] # DNA / RNA

		self.batch_size = 160
		self.max_iter = 50


		if self.hairpin_predictor_name == 'Nupack':
			self.hairpin_predictor = Nupack_Hairpin_Predictor(celsius = self.celsius)
		if self.hairpin_predictor_name == 'Seqfold':
			self.hairpin_predictor = Seqfold_Hairpin_Predictor(celsius = self.celsius)

	def generate(self):
		candidates = np.array([random_seq(self.oligo_length, self.output_format) for i in range (self.batch_size)])
		gc_arr = np.array([GC_content(seq) for seq in candidates])

		gc_cond = np.array(gc_arr > min(self.GC_content_bounds)) & np.array(gc_arr < max(self.GC_content_bounds))
		hairpin_cond = np.array(self.hairpin_predictor.predict(candidates) >= self.hairpin_thr)

		return candidates[gc_cond*hairpin_cond]
		
	def generate_robust(self):
		res = np.empty(shape = 0)
		all_res = np.empty(shape = 0)

		counter = 0
		while len(all_res) < self.num_oligos:
			res = self.generate()
			all_res  = np.unique(np.concatenate((res, all_res), axis = -1))
			counter += 1
			if counter >= self.max_iter:
				OligamaWarning("Max iteration limit reached. Found less oligos than requested", self)
				break

		return all_res[:self.num_oligos]
	
	def describe_oligos(self, seq_arr, toxls = True):
		if len(seq_arr) == 0:
			mes = "No oligos found. Please try less strict conditions. Found less oligos than requested"
			OligamaWarning(mes, self)

		gc_arr = np.array([GC_content(seq) for seq in seq_arr])
		hairpin_arr = np.array(self.hairpin_predictor.predict(seq_arr))
		d_target = {"seq":seq_arr, "GC_content":gc_arr, "Hairpin energy, kJ/mol": hairpin_arr}
		
		df = pd.DataFrame(data = d_target)
		if toxls:
			df_to_excel([df], ['Sheet1'], self.output_folder/"Random_init_table.xlsx")
		return df
