import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from itertools import combinations_with_replacement

from ..data_loaders import Dataloader
from ..utils import df_to_excel


class Affinity_Matrix(Dataloader):
	def __init__(self, input_file, output_folder):
		super().__init__(input_file, output_folder)
		super().add_init()
	def predict(self, units = None):
		if units is None:
			units = self.metric

		n = len(self.target_seqs)
		seq_arr1 = np.repeat(self.target_seqs, n)
		seq_arr2 = np.tile(self.target_seqs, n)
		aff_matr = self.aff_predictor.predict(seq_arr1, seq_arr2, units = units).reshape(n, n)
		
		return aff_matr
	
	def heatmaps_to_png(self, aff_matr, filename_heatmap="Affinity_heatmap.png", filename_hairpin="Hairpin_heatmap.png"):
		plt.figure()
		annot = True
		if len(self.target_seqs) >= 10:
			annot = False
		sns_plot = sns.heatmap(aff_matr,
					mask = np.tril(np.ones_like(aff_matr, dtype=bool), k = -1),
					annot = annot,
					xticklabels = self.target_names,
					yticklabels = self.target_names)
		sns_plot.figure.savefig(self.output_folder/filename_heatmap)

		plt.figure()
		hairpin_en = self.hairpin_predictor.predict(self.target_seqs)
		sns_plot = sns.heatmap(hairpin_en[None, :],
					annot = annot,
					xticklabels = self.target_names,
					yticklabels = False)
		sns_plot.figure.savefig(self.output_folder/filename_hairpin)

	def matrix_to_excel(self, aff_matr, filename="Affinity_matrix.xlsx"):
		hairpin_en = self.hairpin_predictor.predict(self.target_seqs)
		hairpin_df = pd.DataFrame ({"Name": self.target_names, 
									"Seq": self.init_seqs,
									f"Hairpin energy, kJ/mol [{self.hairpin_predictor_name}]": hairpin_en})
		df = pd.DataFrame (aff_matr, columns = self.target_names)
		df.insert(0,'Name','')
		df['Name'] = self.target_names
		df_to_excel([df, hairpin_df], 
					['Affinity', 'Hairpins'],
					self.output_folder/filename)
