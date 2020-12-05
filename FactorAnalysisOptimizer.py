def break_params(theta):
	# TODO
	# MU <- theta()
	# LF <- theta()
	# SIG <- theta()

	MU = None
	LF = None
	SIG = None

	return MU, LF, SIG


def get_f(fam):
	def f(theta):
		MU, LF, SIG = break_params(theta)
		fam.set_theta(MU, LF, SIG)
		return fam.q()
	return f


if __name__ == '__main__':
	# For debugging
	pass
