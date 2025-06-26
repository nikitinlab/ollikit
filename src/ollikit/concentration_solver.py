import numpy as np
import numba as nb

####Функции для нахождения равновесных концентраций
@nb.njit(fastmath=True)
def fpi_step(x, aff_matr_ka, total_conc):
	return total_conc / (1.0 + aff_matr_ka @ x)

@nb.njit(fastmath=True)
def find_equilibrium_conc(total_conc, aff_matr_ka, max_iter = 1000, rtol = 1e-6, atol = 1e-12):
	## Раскомментировать, если вдруг надо работать и с большими диагональными элементами
	# for i in range(aff_matr_ka.shape[0]):
	# 	aff_matr_ka[i, i] *= 2		
	x_new = total_conc.copy()
	x_old = total_conc.copy()
	for k in range(max_iter):
		x_new = fpi_step(x_new, aff_matr_ka, total_conc)
		err = np.abs(x_new - x_old)
		if np.all(err <= atol + rtol*np.maximum(np.abs(x_new), np.abs(x_old))):
			return x_new
		x_old[:] = x_new
	return x_new


##ПРИМЕР использования
# N = 5
# aff_matr_ka = 10**np.random.uniform(4, 12, size = (N, N))
# #Матрица аффинности всегда симметрична
# aff_matr_ka = aff_matr_ka + aff_matr_ka.T
# #На диагонали для простоты пока нули. 
# aff_matr_ka[np.diag_indices(N)] = 0

# total_conc = 10**np.random.uniform(-12, -5, size = N)

# print(aff_matr_ka)
# print(total_conc)
# print(find_equilibrium_conc(total_conc, aff_matr_ka))

	