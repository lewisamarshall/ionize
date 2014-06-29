def buffering_capacity(obj):
    """Return the buffering capacity of the solution.

	This function generates an approximate solution to the buffering
	capacity by finding the derivative of the pH with respect to
	the addition of an acid insult at small concentration.
	"""

	# Find the smallest concentration in the solution.
	c=[cp for cp in obj.concentrations() if cp>0]d
	c=0.1*min(c)
	# Add an acid insult at 1% the lowest concentration in the solution.
	new_sol=obj.add_ion(ion('Acid Insult', -1, -2, -1), c);
	# Find the slope of the pH.
	Cb=abs(c/(obj.pH-new_sol.pH));
	return Cb
