#!/usr/bin/env python
# 12/14/2014
# Charles O. Goddard

import itertools
import sys

from galois import GF, GFpoly, Polynomial


def nth_root_of_unity(GFqn, codeword_len):
	for candidate in map(GFqn, range(2, codeword_len)):
		if (candidate**codeword_len).tonumber() != 1:
			continue
		good = True
		for i in range(1, codeword_len - 1):
			if (candidate**i).tonumber() == 1:
				good = False
				break
		if good:
			return candidate
	raise Exception("Couldn't find a value for alpha!")


if __name__ == '__main__':
	if len(sys.argv) != 4:
		print 'Usage:\n\t%s base power delta\n\n Generate RS code over GF(base**power) with minimum hamming distance delta' % (sys.argv[0], )
		sys.exit(1)

	base = int(sys.argv[1])
	power = int(sys.argv[2])
	delta = int(sys.argv[3])

	if delta > base**power - 1:
		print 'Delta must be <=', base**power - 1
		sys.exit(1)

	# Generate $GF(q^n)$ for $q=base $n$=power and find our root of unity
	GFqn = GFpoly(base, power)
	alpha = nth_root_of_unity(GFqn, base**power - 1)

	# Calculate the generator polynomial $G(x) = \prod_{i=1}^{\delta}{(x - \alpha^i)}$
	Gpoly = Polynomial((GFqn(1),))
	for b in range(1, delta):
		Gpoly *= Polynomial((-(alpha ** b), GFqn(1)))

	codeword_len = base**power - 1
	message_len = base**power - len(Gpoly)

	# Compute the codeword corresponding to the message [1, 2, ... message_len]
	message = Polynomial([GFqn(i + 1) for i in range(message_len)])
	codeword = message * Gpoly

	print 'Alpha =', alpha.tonumber()
	print 'g(x) =', Gpoly
	print 'Message', message.asvector(), 'maps to', codeword.asvector()
	print 'message % g(x) =', message % Gpoly
	print 'codeword % g(x) =', codeword % Gpoly

	maxErrors = (codeword_len - message_len) // 2
	print 'Can detect', codeword_len - message_len, 'errors'
	print 'Can correct', maxErrors, 'errors'
