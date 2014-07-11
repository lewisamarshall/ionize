import numpy

def onsager_fuoss(obj):
    """Return the Onsager-Fuoss correction to the mobilities of ions."""
    # Initialize the empty variables that you will need.
    # Omega = mobility / F / z
    omega = []
    # z_list is an unsorted list of all charges of all ions
    z_list = []
    # conc_list is an unsorted list of all concentrations of all charge states.
    conc_list = []

    # All three share the same order

    # populate them.
    for i in range(len(obj.ions)):
        omega.extend([m/obj._F/z for m, z in
                     zip(obj.ions[i].absolute_mobility, obj.ions[i].z)])
        z_list.extend(obj.ions[i].z)
        conc_list.extend([obj.concentrations[i]*f for f in
                          obj.ions[i].ionization_fraction(obj.pH, obj.I)])

    # add H+ and OH- ions
    omega.extend([obj.H.absolute_mobility/obj.F/1.0,
                  obj.OH.absolute_mobility/obj.F/-1.0])
    z_list.extend([1, -1])
    conc_list.extend([obj.cH, obj.cOH])

    n_states = len(omega)

    omega = numpy.array(omega)
    z_list = numpy.array(z_list)
    conc_list = numpy.array(conc_list)

    # potential is the (chemical?) potential of each ion.
    potential = conc_list*z_list**2/2/obj.I

    # initialize and populate the h matrix.
    h = [[0]*n_states]*n_states
    for j in range(n_states):
        for i in range(n_states):
            try:
                h[i][j] = potential[i]*omega[i]/(omega[i]+omega[j])
            except:
                # Added to avoid divide by zero. Unsure if this is correct.
                h[i][j] = potential[i]/2
    h = numpy.array(h)
    d = numpy.diag(numpy.sum(h, 1))
    B = 2*(h+d)-numpy.identity(n_states)

    r = numpy.zeros(n_states, 6)

    for i in len(r[:, 0]):
        try:
            r[i, 0] = (z_list[i]-sum(z_list*potential) /
                       sum(potential/omega)*(1/omega[i]))
        except ZeroDivisionError:
            r[i, 0] = 0

    for i in range(1, 6):  # used to be 2:6, assume changes based on index.
        r[:, i] = B * r[:, i-1]

    c = [0.2929, -0.3536, 0.0884, -0.0442, 0.0276, -0.0193]
    # coefficients in onsager-fuoss paper

    factor = numpy.dot(c, numpy.transpose(r))

    A_prime = obj.F*0.78420
    B_prime = 31.41e-9

    mob_new = (obj.F*omega-(A_prime*z_list*factor*omega+B_prime)*sqrt(obj.I) /
               (1+1.5*sqrt(obj.I)))

    mob_new = (mob_new*z_list)

    # split the new mobility values into cells that match the molecules
    index = 0
    mobility = [None]*(len(obj.ions)+1)
    for i in range(len(obj.ions)):
        mobility[i] = mob_new[index:(index+length(obj.ions[i].z)-1)]
        index+=length(obj.ions[i].z)

    mobility[-1]=mob_new[index:end]

    return (mobility, omega, z_list, conc_list)
