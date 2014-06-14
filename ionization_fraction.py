function i_frac=ionization_fraction(obj, pH, I ,index)
	% This function takes the ion and the pH and the ionic strength. 
	% It computes the fraction of ion in each of the ionization z states
	% in z. If index is specified, it will only return the ionization 
	% in ionization state z(index). 
			
	% If ionic strength is not specified, set it to zero. 
	if ~exist('I', 'var')
		I=0;
	end
			
	% Sanitize the pH input. 
	if ~isnumeric(pH)
		error('pH should be a number.')
	end
			
	if length(pH)~=1
		pH=pH(1); 
		warning('Ionization fraction only takes a single pH. Using pH(1).')
	end
			
	% Get the vector of products of acidity constants.
    L=obj.L(I);
	% Compute the concentration of H+ from the pH.
    cH=10.^(-pH)/obj.activity_coefficient(I,1);
            
	% Calculate the denominator of the function for ionization fraction.
    i_frac_denom=sum(L.*bsxfun(@power, cH, obj.z0),2);
            
	%Calculate the vector of ionization fractions
    i_frac=L.*bsxfun(@power,  cH, obj.z0)./i_frac_denom;
	i_frac=i_frac(obj.z0~=0);
	% If index is specified, return only the ionization fraction of z(i).
	if exist('index', 'var')
		try
			i_frac=i_frac(index);
		catch
			% Send a warning if the index is meaningless. 
			% Still return the vector. 
			warning('Specified index is out of bounds.')
		end
	end
end