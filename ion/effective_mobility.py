def effective_mobility(obj, pH, I=0):
	'''Returns the effective mobility of the ion at a given pH and I.'''

	if obj.actual_mobility:
		actual_mobility=obj.actual_mobility
	else:
		actual_mobility=obj.robinson_stokes_mobility(I)

	i_frac=ionization_fraction(obj, pH, I)
	effective_mobility=sum([i*j for i,j in zip(i_frac, actual_mobility)])

	return effective_mobility
