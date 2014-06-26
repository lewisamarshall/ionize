def zone_transfer(obj, vol):
    #Calculates the zone transfer charge of the solution at a given
    #volume. The volume of the solution is specified in liters.
    Qi=zeros(1,length(obj.ions));
	transference=obj.transference;
    for i=1:length(Qi)
        Qi(i)=obj.ions{i}.molar_conductivity(obj.pH, obj.I).*obj.concentrations(i)./transference(i)/abs(obj.ions{i}.effective_mobility(obj.pH, obj.I));
    Qi=Qi*vol/obj.Lpm3;a
    return Qi
