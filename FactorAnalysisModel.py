# Module is entitle to the generation of data given
# MU
# LOADING FACTORS
# SIGMA

class FactorAnalysisModel(object):
	def __init__(self, Y):
		self.Y = Y
		self._CZ_MU = None
		self._CZ_LF = None
		self._CZ_SIG = None
		self._J_SIG = None
		self._CY_SIG = None

	def set_theta(self, MU, LF, SIG):
		# TODO:
		self._CZ_MU = MU
		self._CZ_LF = LF
		self._CZ_SIG = SIG
		self._J_SIG = None
		self._CY_SIG = None

	def get_J_mu(self):
		# TODO:
		# return np.array [p + m]
		pass

	def get_J_sig(self):
		# TODO:
		# return np.array [(p + m) x (p + m)]
		if self._J_SIG is None:
			# calculate
			sig = None
			self._J_SIG = sig

		return self._J_SIG

	def get_CY_mu(self):
		# TODO:
		# return np.array [p x 1]
		pass

	def get_CY_sig(self):
		# TODO:
		# return np.array [p x p]
		if self._CY_SIG is None:
			sig = None
			self._CY_SIG = sig

		return self._CY_SIG

	def q_y(self, Y):
		# TODO:
		# Y [m x 1]
		result = None
		return result

	def q(self):
		result = 0
		for y in self.Y:
			result += self.q_y(y)

		return result

if __name__ == '__main__':
	# TODO:
	# For debugging
	pass
