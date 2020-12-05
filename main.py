from FactorAnalysisDataGeneratingProcess import gen_data
from FactorAnalysisOptimizer import get_f
from FactorAnalysisModel import FactorAnalysisModel as FAM

if __name__ == '__main__':
	# MU [m x 1]
	# LF LOADING FACTORS [m x p]
	# SIGMA [mx1] representation of diagonal matrix -> to be turned into proper diag matrix
	# n: number of generated Ys
	# Y ~ N(MU + LF*Z, SIG) [m x n]

	MU = None
	LF = None
	SIGMA = None
	n = None

	Y = gen_data(MU, LF, SIGMA, n)
	fam = FAM(Y)
	f = get_f(fam)

	# TODO: Optimize on f
	pass
