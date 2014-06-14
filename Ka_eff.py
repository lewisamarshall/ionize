function Ka_eff=Ka_eff(obj, I)
	% Uses the ionic strength correction function from
	% Dubye-Huckle theory to calculate the activity coefficients, 
	% and from this, compute the effective Ka values for the ion. 
			
	% If the ionic strength is zero, simply return the Ka's. 
	if I==0;
		Ka_eff=obj.Ka; 
		return
	end
			
	% Make the effective Ka vector the same size as the Ka vector.
	Ka_eff=obj.Ka;
			
	gam_i=obj.activity_coefficient(I);
	gam_h=obj.activity_coefficient(I, 1);
			
	% For each acidity coefficient, get the effective 
	% coefficienty by multiplying by activities.
	for i=1:length(Ka_eff)
		Ka_eff(i)=obj.Ka(i)*gam_i(i+1)/gam_i(i)/gam_h;
	end

end