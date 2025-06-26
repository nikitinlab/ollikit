import scipy.optimize
import numpy as np
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

def create_mat_balance_matrix(n):
    arr = np.array(np.triu_indices(n, 1))
    a = (arr[0][:, np.newaxis] == np.arange(n))
    b = (arr[1][:, np.newaxis] == np.arange(n))
    A_double = (a|b).T.astype(np.int32)

    A_single = np.eye(n)
    return np.concatenate((A_single, A_double), axis = 1)


def find_eq_conc(g, total_conc, celsius):

    n = len(total_conc)
    beta = 1/((273+celsius)*8.31)
    
    g = g * beta
    A = create_mat_balance_matrix(n)

    init_point = np.ones((n**2 + n)//2)*1e-9
    
    cons = scipy.optimize.LinearConstraint(A, lb = total_conc, ub = total_conc)

    hard_bounds = [(0, np.inf)]*len(init_point)

    res = scipy.optimize.minimize(lambda x: np.sum(x*(np.log(x) + g - 1)),
                        init_point,
                        constraints=cons,
                        bounds = hard_bounds,
                        tol = 1e-15, 
                        jac = lambda x: g + np.log(x),
                        )
    return res.x

if __name__ == "__main__":
    n = 3
    g = np.array([0, 0, 0, -52802.10, -15965.45, -42710.79])
    total_conc = np.ones(n)*1e-7
    res = find_eq_conc(g, total_conc)
    print(res)