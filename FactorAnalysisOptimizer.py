import numpy as np

def break_params(theta, y_dim, z_dim):
	# TODO
	# MU <- theta()
	# LF <- theta()
	# SIG <- theta()

	MU = theta[:,y_dim]
	LF = np.matrix(theta[y_dim:-y_dim], (y_dim, z_dim))
	SIG = theta[-y_dim:]

	return MU, LF, SIG


def get_f(fam, y_dim, z_dim):

	def f(theta):
		MU, LF, SIG = break_params(theta, y_dim, z_dim)
		fam.mvd_new.set_theta(MU, LF, SIG)
		return fam.q()

	return f

def get_g(fam, y_dim, z_dim):

	def g(theta):
		MU, LF, SIG = break_params(theta, y_dim, z_dim)
		fam.mvd_old.set_theta(MU, LF, SIG)

		f = get_f(fam, y_dim, z_dim)

		# optimize on f

		return fam.q()

	return g


if __name__ == '__main__':
	# For debugging
	pass
