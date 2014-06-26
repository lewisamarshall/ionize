def equil_offset(obj,I_i):
	''' EQUIL_OFFSET finds the error in the ionic strength.
		Takes an ionic strength, then uses it to calculate
		a new ionic strenght using the new equilibrum coefficents.
		FIND_EQUILIBRIUM finds the root of this function.'''
	pH=obj.calc_pH(I_i);
	I_f=obj.calc_I(pH, I_i);

	Residual=(I_f-I_i);
	return Residual
