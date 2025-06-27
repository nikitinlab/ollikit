import numpy as np
import pandas as pd
import time

from ..data_loaders import Olig_Finder_Dataloader
from ..utils import df_to_excel, multiple_crossover
from ..exceptions import OligamaWarning


class Olig_Finder(Olig_Finder_Dataloader):
	def __init__(self, input_file, output_folder):
		super().__init__(input_file, output_folder)
		self.batch_size = 160

		target_hairpin_en = np.array(self.hairpin_predictor.predict(self.target_seqs))
		if min(target_hairpin_en) < self.hairpin_en_thr:
			mes = f"""Hairpin energy threshold ({self.hairpin_en_thr:.2f} kJ/mol) is greater than hairpin energy for 
			{self.target_names[np.argmin(target_hairpin_en)]} ({(min(target_hairpin_en)):.2f} kJ/mol). Algorithm may not converge."""
			OligamaWarning(mes, self)

	def find_olig(self):	
		seed_targets = self.target_seqs[self.target_aff_high > 0.5]

		if len(seed_targets) == 0:
			seed_targets = self.target_seqs
		candidates = np.array([multiple_crossover(seed_targets) for i in range(self.batch_size)])

		target_aff_arr = self.batch_affinity_predict(self.target_seqs, candidates)
		target_low_cond = np.prod((target_aff_arr - self.target_aff_low[:, None]) > 0, axis = 0,  dtype = bool)
		target_high_cond = np.prod((target_aff_arr - self.target_aff_high[:, None]) < 0, axis = 0,  dtype = bool)

		hairpin_cond = ~np.array(
			self.hairpin_predictor.predict(candidates), dtype = bool)

		res_condition = target_low_cond * target_high_cond * hairpin_cond

		return candidates[res_condition]

	def find_olig_robust(self):
		res = np.empty(shape = 0)
		all_res = np.empty(shape = 0)

		start_time = time.time()
		
		while len(all_res) < self.num_oligos:
			if time.time() - start_time > self.timeout:
				mes = f"Timeout of {self.timeout} seconds reached. Found less oligos than requested. Please try less strict conditions"
				OligamaWarning(mes, self)
				break
				
			res = self.find_olig()
			all_res = np.unique(np.concatenate((res, all_res), axis = -1))

		return all_res[:self.num_oligos]
	
	def describe_oligos(self, seq_arr, toxls = False):
		if len(seq_arr) == 0:
			mes = "No oligos found. Please try less strict conditions"
			OligamaWarning(mes, self)
			return pd.DataFrame()
		target_aff_arr = self.batch_affinity_predict(self.target_seqs, seq_arr)

		d_target = {f"Affinity with {self.target_names[i]}": target_aff_arr[i, :]
									 for i in range(target_aff_arr.shape[0])}
		
		df = pd.DataFrame(data = d_target)
		df.insert(loc=0, column='seq', value=seq_arr)
		df['hairpin'] = self.hairpin_predictor.predict(seq_arr)
		inp_df = pd.DataFrame({"name":self.target_names, "seq":self.init_seqs})
		df = pd.concat((inp_df, df))
		if toxls:
			df_to_excel([df],
					['Sheet1'],
					self.output_folder/"Olig_finder_output.xlsx",
					scientific_format_flag = (self.metric == 'Kd'))
		return df

	def batch_affinity_predict(self, seq1_arr, seq2_arr):
		master_seq_arr = self.arr_product(seq1_arr, seq2_arr)
		aff_arr = self.aff_predictor.predict(master_seq_arr[:, 0], master_seq_arr[:, 1], units = self.metric)
		aff_arr = aff_arr.reshape((len(seq1_arr), len(seq2_arr)))

		return aff_arr

	def arr_product(self, arr1, arr2):
		len1 = len(arr1)
		len2 = len(arr2)
		arr1 = np.array([arr1]*len2).reshape((len2, len1) ).T
		arr2 = np.array([arr2]*len1).reshape((len1, len2) )
		return np.concatenate((arr1[:,:, None], arr2[:,:, None]), axis = 2).reshape(-1, 2)

