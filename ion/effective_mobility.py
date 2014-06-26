def effective_mobility(obj, pH, I=0):
    # This function takes the ion and a pH and an ionic strength. It calls the ionization
    # fraction function and uses this information to compute the
    # effective mobility of the ion.

	if ~isempty(obj.actual_mobility):
		actual_mobility=obj.actual_mobility
	else:
		actual_mobility=obj.robinson_stokes_mobility(I)

	effective_mobility=zeros(size(pH))

    for i = 1:length(pH):
        i_frac=ionization_fraction(obj, pH(i), I)
        effective_mobility(i)=sum(i_frac.*actual_mobility)

    return effective_mobility
