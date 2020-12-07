from FactorAnalysisOptimizer import break_params
import numpy as np

# Module is entitle to the generation of data given
# MU
# LOADING FACTORS
# SIGMA

def gen_data(MU, LF, SIG, n):
	# TODO
	# Y ~ N(MU + LF*Z, SIG)
	# MU
	# LOADING FACTORS
	# SIGMA
	# n: number of generated Ys
	# return a [m x n] matrix
	m = MU.shape[0]
	p = MU.shape[1]

	epsilons = np.random.normal(0, 1, (m, n))
	epsilons = epsilons*SIG.reshape((m, 1))

	Zs = np.random.normal(0, 1, (p, n))

	Ys = MU + LF@Zs + epsilons

	return np.array(Ys)


if __name__ == '__main__':
	# For debugging
	# We define
	# MU
	# LOADING FACTORS
	# SIGMA
	# n: number of generated Ys
	# and cal gen_data function
	# For debugging
	y_dim = 2
	z_dim = 1
	CZ_MU = np.zeros((2, 1))
	CZ_L = np.matrix([[1], [3]])
	CZ_SIG = np.array([1, 1])

	Y = gen_data(CZ_MU, CZ_L, CZ_SIG, 3)
