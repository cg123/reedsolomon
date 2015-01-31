#!/usr/bin/env python
# 12/14/2014
# Charles O. Goddard

class Polynomial(tuple):
	def __new__(cls, iterable):
		if isinstance(iterable, int) and iterable == 0:
			return super(Polynomial, cls).__new__(cls, [])
		t = tuple(iterable)
		while t and t[-1] == 0:
			t = t[:-1]
		return super(Polynomial, cls).__new__(cls, t)

	def __add__(self, rhs):
		if not self:
			return rhs
		res = [type(self[0])(0)]*max(len(self), len(rhs))
		for i in range(len(res)):
			if i < len(self):
				res[i] += self[i]
			if i < len(rhs):
				res[i] += rhs[i]
		return type(self)(res)

	def __sub__(self, rhs):
		if not self:
			return -rhs
		res = [type(self[0])(0)]*max(len(self), len(rhs))
		for i in range(len(res)):
			if i < len(self):
				res[i] += self[i]
			if i < len(rhs):
				res[i] -= rhs[i]
		return type(self)(res)

	def __mul__(self, rhs):
		if not self:
			return self
		res = [type(self[0])(0)]*(len(self) + len(rhs) - 1)
		for i in range(len(rhs)):
			for j in range(len(self)):
				res[i + j] += self[j] * rhs[i]
		return type(self)(res)

	def __neg__(self):
		res = []
		for x in self:
			res.append(-x)
		return type(self)([-x for x in self])

	def divmod(self, rhs):
		if len(self) < len(rhs):
			return type(self)((rhs[0] - rhs[0],)), self

		shiftlen = len(self) - len(rhs)
		N = list(self)
		D = [self[0] - self[0]] * shiftlen + list(rhs)
		
		quot = []
		divisor = D[-1]
		for i in range(shiftlen + 1):
			mult = N[-1] / divisor
			quot = [mult] + quot

			if mult != 0:
			    d = [mult * u for u in D]
			    N = [u - v for u, v in zip(N, d)]

			del N[-1]
			del D[0]
		return type(self)(quot), type(self)(N)

	def __div__(self, rhs):
		return self.divmod(rhs)[0]
	def __mod__(self, rhs):
		return self.divmod(rhs)[1]


	def __pow__(self, x):
		if x == 0:
			return type(self)((type(self[0])(1),))
		res = self
		for i in range(x - 1):
			res = type(self)(res * self)
		return type(self)(res)

	def asvector(self):
		return '[' + ','.join(map(str, self)) + ']'

	def __str__(self):
		if not self:
			return '0'
		res = []
		for i in range(len(self) - 1, -1, -1):
			if self[i] == 0:
				continue
			prefix = ''
			if self[i] > 1:
				prefix = str(self[i]) + ' '

			if i == 0:
				res.append(str(self[i]))
			elif i == 1:
				res.append(prefix + 'x')
			else:
				res.append(prefix + 'x^' +str(i))
		return ' + '.join(res)

	def __repr__(self):
		return self.__str__()

def GF(q, basetype=int):
	class _GFq(basetype):
		def __new__(cls, value):
			while value < 0:
				value += q
			return super(_GFq, cls).__new__(cls, value % q)

		def __add__(self, rhs):
			return _GFq((basetype(self) + rhs))
		def __radd__(self, lhs):
			return _GFq((basetype(self) + lhs))
		def __sub__(self, rhs):
			return _GFq((basetype(self) - rhs))
		def __rsub__(self, lhs):
			return _GFq((lhs - basetype(self)))
		def __mul__(self, rhs):
			return _GFq((basetype(self) * rhs))
		def __rmul__(self, lhs):
			return _GFq((basetype(self) * lhs))
		def __pow__(self, rhs):
			return _GFq((basetype(self) ** rhs))
		def __eq__(self, rhs):
			return basetype(self) == basetype(rhs)

		def __abs__(self):
			return self

		def __div__(self, rhs):
			if rhs == 0:
				raise Exception('Divide by zero')
			remainder = self
			res = _GFq(0)
			while remainder != 0:
				remainder = remainder - rhs
				res += 1
			return res

		def __neg__(self):
			return _GFq(q - int(self))

	return _GFq

def GFpoly(q, n, primpoly=None, basetype=int):
	_GFq = GF(q, basetype)

	def poly_from_idx(pidx):
		p = []
		for _ in range(n + 2):
			p.append(_GFq(pidx % q))
			pidx = pidx // q
		return Polynomial(p)
	def idx_from_poly(p):
		res = 0
		for i in range(len(p)):
			res += int(p[i]) * (q**i)
		return res

	def is_primitive(p, check_irreducible=True):
		if len(p) == 1 or p[0] == 0:
			return False

		if check_irreducible:
			# Check if p is irreducible
			for pidx in range(q, q**n):
				if pidx % q == 0:  # skip multiples of x
					continue
				if not p % poly_from_idx(pidx):
					return False

		for i in range(n, q**n - 1):
			divisor = [0] * (i + 1)
			divisor[0] = _GFq(-1)
			divisor[i] = _GFq(1)
			divisor = Polynomial(divisor)
			if not divisor % p:
				return False
		return True


	def irreducible_sieve():
		is_irreducible = [1] * (q**(n + 1) + 1)
		next_idx = 2
		while next_idx < len(is_irreducible):
			n0 = poly_from_idx(next_idx)
			for mul in range(2, q**n):
				compound_idx = idx_from_poly(n0 * poly_from_idx(mul))
				if compound_idx >= len(is_irreducible):
					continue
				is_irreducible[compound_idx] = 0

			next_idx += 1
			while next_idx < len(is_irreducible) and not is_irreducible[next_idx]:
				next_idx += 1

		for i in range(len(is_irreducible)):
			if is_irreducible[i]:
				yield poly_from_idx(i)

	def find_primitive():
		GFq = GF(q)
		irreducibles = list(irreducible_sieve())
		for i, p in enumerate(irreducibles):
			if len(p) != n + 1:
				continue
			if is_primitive(p):
				return p

		raise Exception("no primitive polynomials found")


	if not primpoly:
		primpoly = find_primitive()
	primpoly = Polynomial(map(_GFq, primpoly))
	print 'primitive polynomial is', primpoly
	class _GFqn(Polynomial):
		def __new__(cls, x):
			if isinstance(x, int):
				res = []
				for i in range(n):
					res.append(_GFq(x % q))
					x = x // q
				return _GFqn(res)

			poly = Polynomial(map(_GFq, x))
			if len(poly) < len(primpoly):
				return super(_GFqn, cls).__new__(cls, poly)
			return super(_GFqn, cls).__new__(cls, poly % primpoly)

		def tonumber(self):
			res = 0
			for i in range(len(self)):
				res += int(self[i]) * (q**i)
			return res

		def __str__(self):
			return str(self.tonumber())

		def __eq__(self, rhs):
			if isinstance(rhs, int):
				return self.tonumber() == rhs
			return super(_GFqn, self).__eq__(rhs)


	return _GFqn
