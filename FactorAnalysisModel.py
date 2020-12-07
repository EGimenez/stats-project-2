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
		self._J_sig_inv = None
		self._J_sig_log_det = None

	def set_theta(self, MU, LF, SIG):
		self._CZ_MU = MU
		self._CZ_LF = LF
		self._CZ_SIG = SIG
		self._J_SIG = None
		self._CY_SIG = None
		self._J_sig_inv = None
		self._J_sig_log_det = None

	def get_X(self):
		n = self.Y.shape[1]
		p = self.get_p()

		z = np.zeros((p, n))

		return np.concatenate((z, self.Y))

	def get_J_mu(self):
		p = self.get_p()
		z = np.zeros((p, 1))
		return np.concatenate((z, self._CZ_MU), axis=0)

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
									sig[i, j] += self._CZ_LF[i-p, k]**2 # ?????
								sig[i, j] += self._CZ_SIG[i-p]**2
							else:
								for k in range(p):
									sig[i, j] += self._CZ_LF[i-p, k]*self._CZ_LF[j-p, k] # ?????

			self._J_SIG = sig

		return self._J_SIG

	def get_CY_mu(self, y):
		p = self.get_p()
		m = self.get_m()
		y = y.reshape((m, 1))
		CZ_MU = self._CZ_MU
		J_SIG = self.get_J_sig()
		sig_12 = J_SIG[:p, p:]
		sig_22 = J_SIG[p:, p:]
		sig_22_inv = np.linalg.inv(sig_22)

		return sig_12@sig_22_inv@(y - CZ_MU)

	def get_J_sig_inv(self):
		if self._J_sig_inv is None:
			self._J_sig_inv = np.linalg.inv(self.get_J_sig())

		return self._J_sig_inv

	def get_CY_sig(self):
		# TODO:
		# return np.array [p x p]
		p = self.get_p()
		m = self.get_m()

		if self._CY_SIG is None:

			J_sig = self.get_J_sig()
			sig11 = J_sig[:p, :p]
			sig12 = J_sig[:p, p:]
			sig21 = J_sig[p:, :p]
			sig22 = J_sig[p:, p:]

			sig = sig11 - sig12 @ np.linalg.inv(sig22) @ sig21
			self._CY_SIG = sig

		return self._CY_SIG

	def get_J_sig_log_det(self):
		if self._J_sig_log_det is None:
			self._J_sig_log_det = np.log(np.linalg.det(self.get_J_sig()))

		return self._J_sig_log_det

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
		# y [m x 1]
		# x [(p + m) x 1]

		result = 0
		p = self.mvd_old.get_p()

		CY_mu_old = self.mvd_old.get_CY_mu(x[p:])
		CY_sig_old = self.mvd_old.get_CY_sig()
		J_mu_new = self.mvd_new.get_J_mu()
		J_sig_new_inv = self.mvd_new.get_J_sig_inv()

		# -1/2 (x-J_mu_new)' * J_sig_new_inv * (x -J_mu_new) -1/2 log(det(J_sig_new))
		# -1/2 log(det(J_sig_new))
		result += -1/2 * self.mvd_new.get_J_sig_log_det()

		# -1/2 (x-J_mu_new)' * J_sig_new_inv * (x -J_mu_new)
		for i in range(J_sig_new_inv.shape[0]):
			for j in range(J_sig_new_inv.shape[1]):
				if (i < p) and (j < p):
					# We are in Z,Z
					result += J_sig_new_inv[i, i] * (CY_sig_old[i, j] + CY_mu_old[i] * CY_mu_old[j])
				elif (i >= p) and (j >= p):
					result += J_sig_new_inv[i, j] * (x[i] - J_mu_new[i])*(x[j] - J_mu_new[j])
				else:
					if i < j:
						result += J_sig_new_inv[i, j] * (x[j] - J_mu_new[j]) * CY_mu_old[i]
					else:
						result += J_sig_new_inv[i, j] * (x[i] - J_mu_new[i]) * CY_mu_old[j]

		return result

	def q(self):
		result = 0
		X = self.mvd_old.get_X()
		for i in range(X.shape[1]):
			result += self.q_x(X[:, i])

		return result

if __name__ == '__main__':
	# For debugging
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

	fam.q()

	print('Hey ho lets go')