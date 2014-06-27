def ionization_fraction(obj, pH, I=0,index):
	'''Return the ionization fractions of an ion at a given pH and
	ionic strength.'''

	# Get the vector of products of acidity constants.
    L=obj.L(I)
	# Compute the concentration of H+ from the pH.
    cH=10.^(-pH)/obj.activity_coefficient(I,1)

	# Calculate the denominator of the function for ionization fraction.
    i_frac_denom=sum(L.*bsxfun(@power, cH, obj.z0),2)

	#Calculate the vector of ionization fractions
    i_frac=L.*bsxfun(@power,  cH, obj.z0)./i_frac_denom
	i_frac=i_frac(obj.z0~=0)
	# If index is specified, return only the ionization fraction of z(i).
	if exist('index', 'var'):
		try
			i_frac=i_frac(index)
		catch
			# Send a warning if the index is meaningless.
			# Still return the vector.
			warning('Specified index is out of bounds.')'''
	return None
	return i_frac
