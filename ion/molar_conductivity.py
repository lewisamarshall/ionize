def molar_conductivity(obj, pH, I=0):
    # This function takes the ion and the pH and computes the molar
    # conductivity. This function can take an ionic strength.
	# Provides conducitivity in Siemens per meter per mole.

	if ~isempty(obj.actual_mobility)
		actual_mobility=obj.actual_mobility;
	else
		actual_mobility=obj.robinson_stokes_mobility(I);
	end

    i_frac=ionization_fraction(obj, pH, I);

    m_conductivity=sum(obj.F*obj.z.*i_frac.*actual_mobility*obj.Lpm3);
    return m_conductivity
