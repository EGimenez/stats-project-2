import numpy as np

# Module is entitle to the generation of data given
# MU
# LOADING FACTORS
# SIGMA

class MultivariateDistribution(object):
	def __init__(self, Y):
		self.Y = Y
		self._CZ_MU = None
		self._CZ_LF = None
		self._CZ_SIG = None
		self._J_SIG = None
		self._CY_SIG = None
		self.J_sig_inv = None

	def set_theta(self, MU, LF, SIG):
		self._CZ_MU = MU
		self._CZ_LF = LF
		self._CZ_SIG = SIG
		self._J_SIG = None
		self._CY_SIG = None
		self.J_sig_inv = None

	def get_X(self):
		n = self.Y.shape[1]
		p = self.get_p()

		z = np.zeros((p, n))

		return np.concatenate((z, self.Y))

	def get_J_mu(self):
		p = self.get_p()
		z = np.zeros((p, 1))
		return np.concatenate(z, self._CZ_MU)

	def get_J_sig(self):
		# return np.array [(p + m) x (p + m)]
		p = self.get_p()
		m = self.get_m()

		if self._J_SIG is None:
			sig = np.zeros((p+m, p+m))

			for i in range(p+m):
				for j in range(p+m):
					if i < p:
						if j < p:
							if i == j:
								sig[i, i] = 1
							else:
								pass
						else:
							sig[i, j] = self._CZ_LF[j-p, i] # ??????
					else:
						if j < p:
							sig[i, j] = self._CZ_LF[i-p, j] # ?????
						else:
							if i == j:
								for k in range(p):
									sig[i, j] += self._CZ_LF[i-p, k]^2 # ?????
								sig[i, j] += self._CZ_SIG[i-p, i-p]
							else:
								for k in range(p):
									sig[i, j] += self._CZ_LF[i-p, k]*self._CZ_LF[j-p, k] # ?????

			self._J_SIG = sig

		return self._J_SIG

	def get_J_sig_inv(self):
		if self.J_sig_inv is None:
			self.J_sig_inv = np.linalg.inv(self.J_sig())

		return self.J_sig_inv

	def get_CY_sig(self):
		# TODO:
		# return np.array [p x p]
		p = self.get_p()
		m = self.get_m()

		if self._CY_SIG is None:

			J_sig = self.get_J_sig()
			sig11 = J_sig[:p, :p]
			sig12 = J_sig[m:, :p]
			sig21 = J_sig[:p, m:]
			sig22 = J_sig[m:, m:]

			sig = sig11 - sig12 * np.linalg.inv(sig22) * sig21
			self._CY_SIG = sig

		return self._CY_SIG

	def get_p(self):
		return self._CZ_LF.shape[1]

	def get_m(self):
		return self._CZ_LF.shape[0]

class FactorAnalysisModel(object):
	def __init__(self, mvd_old, mvd_new):
		self.mvd_old = mvd_old
		self.mvd_new = mvd_new

	def set_theta(self, MU, LF, SIG):
		self.mvd_new.set_theta(MU, LF, SIG)

	def q_x(self, x):
		# TODO:
		# y [m x 1]
		# x [(p + m) x 1]

		result = 0
		p = self.mvd_old.get_p()

		CY_sig_old = self.mvd_old.get_CY_sig()
		J_sig_new_inv = self.mvd_new.get_J_sig_inv()
		J_mu_new = self.mvd_new.get_J_mu()

		# -1/2 (x-mu)' * J_sig_new_inv * (x -mu) -1/2 log(det(J_sig)) * old_distrib
		# -1/2 log(det(J_sig))
		result += self.mvd_new.get_J_sig_log_det()

		# -1/2 (x-mu)' * J_sig_new_inv * (x -mu)
		for i in range(J_sig_new_inv.shape[0]):
			for j in range(J_sig_new_inv.shape[1]):
				if i < p:
					if j < p:
						# We are in Z,Z
						if i == j:
							# We are in Zi, Zi the variance of Zi
							result += CY_sig_old[i, i]*J_sig_new_inv[i, i]
						else:
							# We are in teh Zi, Zj so its zero
							pass
					else:
						# We are in Z,Y the expected value of Z is zero
						pass
				else:
					if j < (p - 1):
						# We are in Y, Z the expected value of Z is zero
						pass
					else:
						# We are in Y, Y
						result += (x[i] - J_mu_new[i])*J_mu_new[i, j]*(x[j] - J_mu_new[j])

		return result

	def q(self):
		result = 0
		X = self.mvd_old.get_X()
		for i in X.shape[1]:
			result += self.q_x(X[:, i])

		return result

if __name__ == '__main__':
	# TODO:
	# For debugging
	pass
