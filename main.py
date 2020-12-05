from FactorAnalysisDataGeneratingProcess import gen_data
from FactorAnalysisOptimizer import get_g
from FactorAnalysisModel import FactorAnalysisModel as FAM
from FactorAnalysisModel import MultivariateDistribution as MVD

if __name__ == '__main__':
	# MU [m x 1]
	# LF LOADING FACTORS [m x p]
	# SIGMA [mx1] representation of diagonal matrix -> to be turned into proper diag matrix
	# n: number of generated Ys
	# Y ~ N(MU + LF*Z, SIG) [m x n]

	CZ_MU = None
	CZ_LF = None
	CZ_SIGMA = None
	n = None

	Y = gen_data(CZ_MU, CZ_LF, CZ_SIGMA, n)
	init_values = None
	mvd_old = MVD(Y, init_values)
	mvd_new = MVD(Y, init_values)
	fam = FAM(mvd_old, mvd_new)
	g = get_g(fam)

	# TODO: Optimize on f
	pass
