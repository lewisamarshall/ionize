function [mobility, omega, z_list, conc_list]=onsager_fuoss(obj)
			
	% Initialize the empty variables that you will need.
	% Omega = mobility / F / z
	omega=[];
	% z_list is an unsorted list of all charges of all ions
	z_list=[];
	% conc_list is an unsorted list of all concentrations of all charge states. 
	conc_list=[];
			
	%All three share the same order
			
	% populate them.
	for i=1:length(obj.ions)
		omega=cat(2, omega, (obj.ions{i}.absolute_mobility./obj.F./obj.ions{i}.z));
		z_list=cat(2, z_list, (obj.ions{i}.z));
		conc_list=cat(2, conc_list, (obj.concentrations(i).*obj.ions{i}.ionization_fraction(obj.pH, obj.I)));	
	end
	
	% add H+ and OH- ions
	omega=cat(2, omega, [obj.H.absolute_mobility, obj.OH.absolute_mobility]./obj.F./[1, -1]); 
	z_list=cat(2, z_list, [1, -1]);
	conc_list=cat(2, conc_list, [obj.cH, obj.cOH]);
	
	n_states=length(omega);
	
	% potential is the (chemical?) potential of each ion.
	potential=conc_list.*z_list.^2./(2*obj.I); 
	% initialize and populate the h matrix. 
    h=zeros(n_states);
	for j=1:n_states
	    for i=1:n_states
		   	h(i,j)=potential(i)*omega(i)/(omega(i)+omega(j));
			% Added to avoid NaN. Unsure if this is the correct strategy.
			if isnan(h(i,j))
				h(i,j)=potential(i)/2;
			end
	    end
	end

	d=diag(sum(h,2));
    B=2*(h+d)-eye(n_states);

	r=zeros(n_states, 6);

    r(:,1) = (z_list-sum(z_list.*potential) / sum(potential.*abs(z_list)./omega./abs(z_list)) * (1./omega))';
	
	% if r==NaN, r=0. 
	for i=1:n_states
		if isnan(r(i, 1))
			r(i,1)=0;
		end
	end
		
    for i=2:6                                        
    	r(:,i)=B*r(:,i-1);
    end

	c=[0.2929 -0.3536 0.0884 -0.0442 0.0276 -0.0193];
	%coefficients in onsager-fuoss paper

    factor=(c*r');

	A_prime=obj.F*0.78420;
	B_prime=31.41e-9;
	
	mob_new=obj.F*omega-(A_prime  *z_list.*factor.*omega+B_prime).*sqrt(obj.I)./(1+1.5*sqrt(obj.I));

	mob_new=(mob_new.*z_list);
	
	% split the new mobility values into cells that match the
	% molecules
	index=1;
	mobility=cell(1, length(obj.ions)+1);
	for i=1:length(obj.ions)
		mobility{i}=mob_new(index:(index+length(obj.ions{i}.z)-1));
		index=index+length(obj.ions{i}.z);
	end
	mobility{end}=mob_new(index:end);

end