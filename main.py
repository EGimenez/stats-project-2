from FactorAnalysisDataGeneratingProcess import gen_data
from FactorAnalysisModel import FactorAnalysisModel as FAM
from FactorAnalysisModel import MultivariateDistribution
from FactorAnalysisOptimizer import get_f, break_params
import numpy as np
from scipy.optimize import minimize

if __name__ == '__main__':
	# MU [m x 1]
	# LF LOADING FACTORS [m x p]
	# SIGMA [mx1] representation of diagonal matrix -> to be turned into proper diag matrix
	# n: number of generated Ys
	# Y ~ N(MU + LF*Z, SIG) [m x n]

	y_dim = 2
	z_dim = 1
	CZ_MU = np.zeros((2, 1))
	CZ_L = np.matrix([[1], [3]])
	CZ_SIG = np.array([1, 1])

	Y = gen_data(CZ_MU, CZ_L, CZ_SIG, 10)

	mvd_old = MultivariateDistribution(Y)
	mvd_new = MultivariateDistribution(Y)

	init_guess = np.array([0, 0, 0, 0, .9, .9])
	mu_init, l_init, sig_init = break_params(init_guess, y_dim, z_dim)

	mvd_old.set_theta(mu_init, l_init, sig_init)
	mvd_new.set_theta(mu_init, l_init, sig_init)

	fam = FAM(mvd_old, mvd_new)

	f = get_f(fam, y_dim, z_dim)

	res = minimize(f,
	               init_guess,
	               method='L-BFGS-B',
	               options={'maxiter': 1000})


	for s in range(100):
		mu, lf, sig = break_params(res.x, y_dim, z_dim)
		fam.mvd_old.set_theta(mu, lf, sig)
		f = get_f(fam, y_dim, z_dim)
		res = minimize(f,
		               res.x,
		               method='L-BFGS-B',
		               options={'maxiter': 1000})

	mu, lf, sig = break_params(res.x, y_dim, z_dim)

	print('Hey ho lets go')