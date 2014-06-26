def Kw_eff(obj, I=obj.I):
	# Gets the effective water dissociation constant
	# due to activity corrections to H+ and OH-.
	if ~exist('I', 'var')
		I=obj.I;

	gam_h=obj.H.activity_coefficient(I,1);
	Kw_eff=obj.Kw/gam_h^2;
	return Kw_eff
