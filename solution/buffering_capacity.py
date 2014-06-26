def buffering_capacity(obj):
	'''BUFFERING_CAPACITY finds the buffering capacity of a solution.
		This function generates an exact solution to the buffering
		capacity by finding the derivative of the pH with respect to
		the addition of an acid.'''

	# Find the smallest concentration in the solution.
	c=obj.concentrations(obj.concentrations>0);
	c=0.1*min(c);
	# Add an acid insult at 1% the lowest concentration in the solution.
	new_sol=obj.add_ion(ion('Acid Insult', -1, -2, -1), c);
	# Find the slope of the pH.
	Cb=abs(c/(obj.pH-new_sol.pH));
	return Cb
