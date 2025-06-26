import nupack
import numpy as np
from seqfold import dg
from joblib import Parallel, delayed

class CNN_Hairpin_Predictor():
	def __init__(self, celsius = 25):
		self.celsius = celsius
	def predict(self, seq_arr):
		return np.zeros(len(seq_arr))
	


class Nupack_Hairpin_Predictor():
	def __init__(self, celsius = 25):
		self.celsius = celsius
		self.model = nupack.Model(material = 'dna', celsius = celsius)

	def check_hairpin(self, seq):
		s1 = nupack.Strand(seq, name='s1')
		complexes = nupack.SetSpec(max_size=2)
		complex_set = nupack.ComplexSet(strands=[s1], complexes=complexes)
		result = nupack.complex_analysis(complexes=complex_set, model=self.model, compute=['pfunc', 'mfe'])
		energies = [s1_structure.energy for s1_structure in result['(s1)'].mfe]
		return min(energies)

	def predict(self, seq_arr, n_jobs = -1):
		return np.array(
			Parallel(n_jobs = n_jobs) (delayed(self.check_hairpin) (seq) for seq in seq_arr)
		)





class Seqfold_Hairpin_Predictor():
	def __init__(self, celsius = 25):
		self.celsius = celsius


	def predict(self, seq_arr, n_jobs = -1):
		res =  np.array(
			Parallel(n_jobs = n_jobs) (delayed(dg) (seq, self.celsius) for seq in seq_arr)
		)
		res[res > 0] = 0
		return res



