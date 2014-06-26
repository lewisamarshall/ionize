def calc_I(obj, pH, I_guess=0):
    '''CALC_I calculates the ionic strength of a solution.
        CALC_I should only be used during solution initialization,
        after that, the equilibrium ionic strength is stored in obj.I.
        pH must be supplied. If a guess for ionic strength is not supplied,
        ionic strength corrections will not be used.'''


	# Set the ionic strength to zero to start with. It will be counted for each ion.
	I=0;

	# For each ion, add the contribution to ionic strength to the sum.
	for i=1:length(obj.ions);
		I=I+obj.concentrations(i)*sum(((obj.ions{i}.z).^2).*obj.ions{i}.ionization_fraction(pH, I_guess));
	end

	# Add the ionic strength due to water dissociation.
	I=I+obj.cH(pH) + obj.cOH(pH, I_guess);
	# Divide by 2 to get the correct value.
	I=I/2;
    return I
