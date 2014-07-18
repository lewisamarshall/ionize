def zone_transfer(obj, vol):
    """Return the zone transfer charge of the solution at a given volume.

    The volume of the solution is specified in liters.
    """
    Qi = [0]*len(obj.ions)
    transference = obj.transference()
    for i in range(len(Qi)):
        Qi[i] = (obj.ions[i].molar_conductivity(obj.pH, obj.I) *
                 obj.concentrations[i]/transference[i] /
                 abs(obj.ions[i].effective_mobility(obj.pH, obj.I)))
    Qi = [Qp*vol/obj._Lpm3 for Qp in Qi]
    return Qi
