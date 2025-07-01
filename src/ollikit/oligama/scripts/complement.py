from pathlib import Path

import pandas as pd
from ..utils import df_to_excel, compl
from ..data_loaders import Complement_Dataloader

class Complement(Complement_Dataloader):
	def __init__(self, input_data, output_folder):
		super().__init__(input_data, output_folder)

	def predict(self):
		"""
		Генерирует комплементарные и обратные последовательности.

		Returns:
			tuple: (rev_arr, rev_compl_arr, dir_compl_arr) — обратные, обратные комплементарные, прямые комплементарные последовательности.
		"""
		# Генерируем комплементарные последовательности с учетом материала
		rev_compl_arr = [compl(seq, material=material.lower()) 
						for seq, material in zip(self.target_seqs, self.sequence_types)]
		dir_compl_arr = [seq[::-1] for seq in rev_compl_arr]
		rev_arr = [seq[::-1] for seq in self.target_seqs]
		return rev_arr, rev_compl_arr, dir_compl_arr

	def save_to_excel(self, rev_arr, rev_compl_arr, dir_compl_arr, filename="Complement_table_output.xlsx"):
		df = pd.DataFrame({"Name":self.target_names,
							"Initial sequence (5'-3')" : self.target_seqs, 
							"Reverse sequence (3'-5')" : rev_arr,
							"Reverse complement (5'-3')": rev_compl_arr,
							"Complement (3'-5')": dir_compl_arr})
		df_to_excel([df], ['Sheet1'], self.output_folder / Path(filename))

