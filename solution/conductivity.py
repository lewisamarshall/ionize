function conductivity=conductivity(obj)
    % CONDUCTIVITY calcualtes the electrical conductivity of the solution. 
	%	Relies on molar conductivity calculations from ion and total conductivity 
	%	of H+ and OH-. 
			
    conductivity=0;
    for i=1:length(obj.concentrations)
        conductivity=conductivity+...
			obj.ions{i}.molar_conductivity(obj.pH, obj.I) * obj.concentrations(i);
    end
			
    conductivity=conductivity + obj.OH_conductivity; 	% Add contribution for hydroxyl conductivity
    conductivity=conductivity + obj.H_conductivity;		% Add contribution for hydronium conductivity
end