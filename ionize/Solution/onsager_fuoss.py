import numpy
from math import sqrt


def onsager_fuoss(self):
    """Return the Onsager-Fuoss correction to the mobilities of ions.

    This function returns a list of all corrected actual mobilities.
    These mobilities are automatically assigned to the correct ions when
    a solution is initialized."""
    # Initialize the empty variables that you will need.
    # Omega = mobility / F / z
    omega = []
    # z_list is an unsorted list of all charges of all ions
    z_list = []
    # conc_list is an unsorted list of all concentrations of all charge states.
    conc_list = []

    # All three share the same order

    # populate them.
    for i in range(len(self.ions)):
        omega.extend([m/self._F/z for m, z in
                     zip(self.ions[i].absolute_mobility, self.ions[i].z)])
        z_list.extend(self.ions[i].z)
        conc_list.extend([self.concentrations[i]*f for f in
                          self.ions[i].ionization_fraction(self.pH, self.I)])

    # add H+ and OH- ions
    omega.extend([self._H.absolute_mobility[0]/self._F/1.0,
                  self._OH.absolute_mobility[0]/self._F/-1.0])
    z_list.extend([1, -1])
    conc_list.extend([self.cH(), self.cOH()])

    # n_states is the total number of ionic species that we are tracking.
    n_states = len(omega)

    # turn all variables into numpy arrays for matrix math
    omega = numpy.array(omega)
    z_list = numpy.array(z_list)
    conc_list = numpy.array(conc_list)

    # potential is the (chemical?) potential of each ion.
    potential = conc_list*z_list**2/2/self.I

    # initialize and populate the h matrix.
    h = [[0]*n_states]*n_states
    for j in range(n_states):
        for i in range(n_states):
            try:
                h[i][j] = potential[i]*omega[i]/(omega[i]+omega[j])
            except:
                # Added to avoid divide by zero. Unsure if this is correct.
                h[i][j] = potential[i]/2
    # turn h into a numpy array
    h = numpy.array(h)
    d = numpy.diag(numpy.sum(h, 1))
    B = 2*(h+d)-numpy.identity(n_states)

    r = numpy.zeros([n_states, 6])

    for i in range(len(r[:, 0])):
        try:
            r[i, 0] = (z_list[i]-sum(z_list*potential) /
                       sum(potential/omega)*(1/omega[i]))
        except ZeroDivisionError:
            r[i, 0] = 0

    for i in range(1, 6):  # used to be 2:6, assume changes based on index.
        r[:, i] = numpy.dot(B, r[:, i-1])

    # coefficients in onsager-fuoss paper
    c = [0.2929, -0.3536, 0.0884, -0.0442, 0.0276, -0.0193]

    factor = numpy.dot(c, numpy.transpose(r))

    T = self.T
    T_ref = 25
    d = self._dielectric(T)
    d_ref = self._dielectric(T_ref)

    # New temperature corrected coefficients.
    # Temperature corrections are based on the reference values.
    A_prime = self._F*0.78420*((T_ref+273.15)*d_ref/(T+273.15)/d)**(-1.5)
    B_prime = 31.410e-9 * ((T_ref+273.15)*d_ref/(T+273.15)/d)**(-0.5) *\
        self._viscosity(T_ref)/self._viscosity(T)

    # A_prime = self._F*0.78420
    # B_prime = 31.41e-9

    mob_new = (self._F*omega-(A_prime*z_list*factor*omega+B_prime)*sqrt(self.I) /
               (1+1.5*sqrt(self.I)))*z_list

    # transfer the new mobilities from a numpy array back to alist.
    mob_new = mob_new.tolist()

    # split the new mobility values into cells that match the molecules
    index = 0
    mobility = [None]*(len(self.ions)+1)
    for i in range(len(self.ions)):
        mobility[i] = mob_new[index:(index+len(self.ions[i].z))]
        index += len(self.ions[i].z)

    mobility[-1] = mob_new[index:]
    return mobility
