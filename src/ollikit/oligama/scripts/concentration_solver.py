import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

from ..scripts.affinity_matrix import Affinity_Matrix
from ..utils import df_to_excel
from ..data_loaders import *
from ..complex_concentration_minimal import find_eq_conc

 
class Concentration_Solver(Concentration_Solver_Dataloader):
	def __init__(self, input_file, output_folder):
		super().__init__(input_file, output_folder)

		Aff_matr = Affinity_Matrix(input_file, output_folder)
		self.aff_matr = Aff_matr.predict(visualize = False, units = 'gibbs')

		n = len(self.target_seqs)
		self.linear_g = np.concatenate((np.zeros(n), self.aff_matr[np.triu_indices(n, 1)]))
  
		self.inp_conc = np.geomspace(self.input_seq_bounds[0], self.input_seq_bounds[1], self.n_points)
		self.out_conc = np.zeros(self.n_points)

		for ind, name in enumerate(self.target_names):
			if name == self.output_seq_name:
				self.output_ind = ind
			if name == self.input_seq_name:
				self.input_ind = ind

	def predict(self, toxls = True):
		for ind, inp_conc in enumerate(self.inp_conc):

			self.total_conc[self.input_ind] = inp_conc

			ccc = find_eq_conc(self.linear_g, self.total_conc, self.celsius)[self.output_ind]

			self.out_conc[ind] = ccc

		fig = plt.figure(figsize = (7, 5))
		plt.plot(self.inp_conc, self.out_conc)
		plt.scatter(self.inp_conc, self.out_conc)

		ax = plt.gca()
		ax.set_xscale('log')
		ax.set_yscale('log')
		ax.set_xlabel(f"{self.input_seq_name} concentration, M")
		ax.set_ylabel(f"{self.output_seq_name} concentration, M")
		plt.savefig(self.output_folder / Path("Concentration_solver_graph.png"))

		inp_df = pd.DataFrame({"name":self.target_names, "seq":self.target_seqs, "concentration, M":self.total_conc})
		
		df = pd.DataFrame({f"{self.input_seq_name} concentration, M": self.inp_conc,
					 	   f"{self.output_seq_name} concentration, M": self.out_conc})
		if toxls :
			df_to_excel([df, inp_df], ['Curve', "Input"],self.output_folder / Path("Concentration_solver_table.xlsx"), scientific_format_flag = True)
		return df, inp_df #self.out_conc
