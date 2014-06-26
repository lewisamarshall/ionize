def transference(obj):
    '''TRANSFERENCE Gets the fraction of charge carried by each of the ions.
    	Should not precisely add to 1, because some charge is carried
    	by protons and hydroxyls.'''
    T=zeros(1,length(obj.ions));
    for i=1:length(T)
        T(i)=obj.ions{i}.molar_conductivity(obj.pH, obj.I).*obj.concentrations(i);

    T=T/obj.conductivity;
    return T
