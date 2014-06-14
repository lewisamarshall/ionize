function gamma=activity_coefficient(obj, I, z)
	if ~exist('z', 'var')
		z=obj.z0;
	end
			
	% There are two coefficients that are used repeatedly.
	% Specified in Bahga. 
	A=obj.Adh*sqrt(I)/(1+obj.aD*sqrt(I));
	B=obj.Adh*0.1*I; %altered to match spresso code, Adh may not belong
			
	% Use them to calculate the activity coefficients. 
	% These coefficients are for z=+-1, for H+ and OH-
	gamma=z.^2*(B-A);
	gamma=10.^gamma;
end