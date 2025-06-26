import numpy as np


def triangle_to_matr(vec):
	n = int((1 + np.sqrt(1 + 8*len(vec)))//2)
	matr = np.zeros((n, n))
	matr[np.triu_indices(n, 1)] = vec
	return matr + matr.T


def detect_unit(arr):
	if isinstance(arr, list):
		arr = np.array(arr)
	if arr.min() < 0:
		return "energy"
	elif arr.max() > 1:
		return "Ka"
	elif arr.min() < 1e-4 and arr.max() < 0.1:
		return "Kd"
	else:
		return "fraction"

def convert(arr, to_unit = 'fraction', from_unit = None,  temp = 25):
	to_unit = to_unit.lower()
	if isinstance(arr, list):
		arr = np.array(arr)
	if isinstance(arr, float) or isinstance(arr, int):
		arr = np.array([arr])
		
	if from_unit is None:
		from_unit = detect_unit(arr)
	from_unit = from_unit.lower()
	if from_unit == "fraction":
		mask = (arr != 0)
		Ka = np.zeros_like(arr)
		Ka[mask] = 1e6*arr[mask]/(1 - arr[mask])**2
	elif from_unit == "kd":
		Ka = 1/arr
	elif from_unit == "energy" or from_unit == "gibbs":
		Ka = np.exp(-arr*1e3/(8.31*(temp + 273)))
	elif from_unit == "ka":
		Ka = arr
	else:
		print("Unknown FROM unit")
		return None
	
	if to_unit == "fraction":
		mask = (Ka > 1e1)
		res = np.zeros_like(Ka)
		res[mask] = (2*Ka[mask] + 1e6 - np.sqrt(4*Ka[mask]*1e6 + 1e12)) / (2*Ka[mask])
		return res
	elif to_unit == "kd":
		Kd = np.ones_like(Ka)
		Kd[Ka != 0] = 1/Ka[Ka != 0]
		return Kd
	elif to_unit == "energy" or to_unit == "gibbs":
		res = np.zeros_like(Ka)
		res[Ka != 0] = -8.31*(temp + 273)*np.log(Ka[Ka != 0])/1000
		return res
	elif to_unit == "ka":
		return Ka
	else:
		print("Unknown TO unit")
		return None
	

	

	
