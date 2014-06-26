def find_equilibrium(obj):
	''' FIND_EQUILIBRIUM finds the equilibrium ionic strength and pH.
		It uses the fzero root finder to find the equilibrium
		pH and ionic strength of a solution, using the ionic-
		strength-adjusted activity coefficients. This function
		is called when the object is initialized.'''

	# Generate an initial ionic strength guess without using activity corrections
	I=obj.calc_I(obj.calc_pH);
	# Iterate to find the true ionic strength.

	OPTION=optimset('TolX', 1e-6);
	[I,fval,exitflag]=fzero(@obj.equil_offset, I, OPTION);

	if exitflag~=1
		error('Could not find equilibrium.')
	end
	# Use this final ionic strength to find the correct pH.
	pH=obj.calc_pH(I);
	return pH, I
