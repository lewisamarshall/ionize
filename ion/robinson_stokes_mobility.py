def robinson_stokes_mobility(obj, I):
	'''If only an ionic strength is specified, use the Robinson-Stokes
	correction to calculate a new fully ionized mobility'''

	# If a solution object is supplied, use the full onsager fouss correction.

	if isnumeric(I) and  I>=0 and isvector(I) and length(I)==1:
		# Currently using the ionic strength where Bahga 2010
		# uses twice the ionic strength. This appears to work, and follows the SPRESSO implimentation.
		# Likely typo in paper.
		A=0.2297;
		B=31.410e-9;
		#actual_mobility=obj.absolute_mobility-(A.*obj.absolute_mobility+B.*sign(obj.z))*sqrt(I)/(1+obj.aD*sqrt(I));
		actual_mobility=0
	else:
		error('Ionic strength must be specified as a scalar positive value.')

	return actual_mobility
