import numpy as np
from FactorAnalysisModel import FactorAnalysisModel, MultivariateDistribution
from scipy.optimize import minimize


def break_params(theta, y_dim, z_dim):
	# TODO
	# MU <- theta()
	# LF <- theta()
	# SIG <- theta()

	MU = theta[:y_dim].reshape((y_dim, 1))
	LF = theta[y_dim:-y_dim].reshape((y_dim, z_dim))
	SIG = theta[-y_dim:]**2 # + 0.00000001 # This **2 is to avoid restrictions on the optimizer
	return MU, LF, SIG


def get_f(fam, y_dim, z_dim):

	def f(theta):
		MU, LF, SIG = break_params(theta, y_dim, z_dim)
		fam.mvd_new.set_theta(MU, LF, SIG)
		return -fam.q()

	return f


if __name__ == '__main__':
	# For debugging
	y_dim = 2
	z_dim = 1
	CZ_MU_old = np.zeros((2,1))
	CZ_L_old = np.matrix([[1], [3]])
	CZ_SIG_old = np.array([1, 1])
	CZ_MU_new = np.zeros((2,1))
	CZ_L_new = np.matrix([[1], [3]])
	CZ_SIG_new = np.array([1, 1])

	Y = np.random.random((2, 10))

	mvd_old = MultivariateDistribution(Y)
	mvd_old.set_theta(CZ_MU_old, CZ_L_old, CZ_SIG_old)
	mvd_new = MultivariateDistribution(Y)
	mvd_new.set_theta(CZ_MU_new, CZ_L_new, CZ_SIG_new)

	fam = FactorAnalysisModel(mvd_old, mvd_new)

	f = get_f(fam, 2, 1)

	res = minimize(f,
	               np.array([0, 0, 0, 0, .9, .9]),
	               method='L-BFGS-B',
	               options={'maxiter': 1000})

	mu, lf, sig = break_params(res.x, 2, 1)

	print('Hey ho lets go')
