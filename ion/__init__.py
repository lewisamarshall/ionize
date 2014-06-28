import warnings
from math import sqrt

class ion:

    """Describe an ion dissolved in aqueous solution.

    This class draws significantly on Bagha Electrophoresis 2010
    "Ionic strength effects on electrophoretic focusing and separations"
    The ion is defined by a name, and a set of charge states (z).
    Each charge state must have an associated acidity constant (pKa).
    Each charge state must have an associated fully ionized mobility
    (absolute_mobility).

    This is a direct port of the Matlab code written by Lewis Marshall.
    """

    # Weakly private variables
    # These are constants and should not change.
    # Eventually, T may be  removed from the constants list.
    _F = 96485.3415     # Faraday's const.[C/mol]
    _Lpm3 = 1000        # Conversion from liters to m^3
    _T = 298            # Temperature, in Kalvin
    # The following are constants in eqtn 6 of Bahga 2010.
    _Adh = 0.5102  	 # L^1/2 / mol^1/2, approximate for RT
    _aD = 1.5  	     # mol^-1/2 mol^-3/2, approximation

    def __init__(self, name, z, pKa, absolute_mobility):
        """Initialize an ion object."""
        self.name = name
        self.z = z
        self.pKa = pKa
        self.absolute_mobility = absolute_mobility  # Expected in m^2/V/s.
        self.actual_mobility = [None] * len(z)   # Intended to be filled by solution


    # Initialize the object, with checks on the form of the variables.
        # if(nargin ==4 ):
        #     # check that the name is a string
        #     if ischar(name)
        #         obj.name=name;
        #     else
        #         error('Ion name must be a string.')
        #     end
        #
        #     % Check that z is a vector of integers
        #     if all(z==int8(z)) && isvector(z)
        #         obj.z=double(z);
        #     else
        #         error('Charge states must be a vector of integers.')
        #     end
        #
        # 	% Check that the pKa is a vector of numbers of the same length as z.
        #     if isvector(pKa) && length(obj.z)==length(pKa) &&isnumeric(absolute_mobility)
        #         obj.pKa=double(pKa);
        #     else
        #         error ('pKas must be a numeric vector the same size as charge vector.')
        #     end
        #
        # 	%Check that the fully ionized mobility is a vector of numbers the same size as z.
        #     if isvector(absolute_mobility) && length(obj.z)==length(absolute_mobility) && isnumeric(absolute_mobility)
        #         obj.absolute_mobility=double(absolute_mobility);
        #     else
        #         error ('Fully ionized mobilities must be a numeric vector the same size as charge vector.')
        #     end
        #
        # 	% Force the sign of the fully ionized mobilities to match the sign of the charge.
        # 	% This command provides a warning, which you can suppress, with, for example,
        # 	% warning('off','all');
        #     if ~all(sign(obj.z)==sign(obj.absolute_mobility))
        #         obj.absolute_mobility=abs(obj.absolute_mobility).*double(sign(obj.z));
        #         warning('Forcing fully ionized mobility signs to match charge signs.')
        #     end
        #
        # else
        # 	% If there are not four inputs, the initialization fails.
        #     error ('Not enough inputs to describe an ion.')
        # end
        #
        # % After storing the ion properties, ensure that the properties are sorted in order of charge.
        # % All other ion methods assume that the states will be sorted by charge.
        self=self.z_sort();

    def z_sort(obj):
        """Sort the charge states from lowest to highest."""
        # Zip the lists together and sort them by z.
        obj.z, obj.pKa, obj.absolute_mobility=zip(*sorted(zip(obj.z, obj.pKa, obj.absolute_mobility)))
        obj.z=list(obj.z)
        obj.pKa=list(obj.pKa)
        obj.absolute_mobility=list(obj.absolute_mobility)


        # This section will check each charge state to see if it is complete.
        # That is, if there is a charge state -2, there must be a charge state
        # -1. if there is a charge state +3, there must be a +2 and a +1.
        # warn=0
        # for i in range(len((obj.z))):
        #     if obj.z(i) < -1 and obj.z(i+1)!=obj.z(i)+1:
        #         warn=1
        #     elif obj.z(i)>1 and obj.z(i-1)!=obj.z(i)-1:
        #         warn=1
        #
        #     #Send a single warning if any charge state is missing
        # if warn:
        #     warnings.warn('Charge states missing.')
        return obj

    def Ka(obj):
        """Return the Kas based on the pKas."""
        Ka=[10.**-p for p in obj.pKa]
        return Ka

    def z0(obj):
        """Return the list of charge states with 0 inserted."""
        z0=[0]+obj.z;
        z0=sorted(z0);
        return z0

    from ionization_fraction import ionization_fraction
    from activity_coefficient import activity_coefficient
    from effective_mobility import effective_mobility
    from Ka_eff import Ka_eff
    from L import L
    from molar_conductivity import molar_conductivity
    from robinson_stokes_mobility import robinson_stokes_mobility

if __name__=='__main__':
    hcl=ion('hydrochloric acid', [-1, -2], [-2, 3], [76, 89])
    print hcl
    print hcl.name
    print hcl.z
    print hcl.pKa
    print hcl.absolute_mobility
    print hcl.robinson_stokes_mobility(.1)
    print hcl.Ka()
    print hcl.z0()
    #  print hcl.ionization_fraction(7)
    print hcl.activity_coefficient(.03)
    help(ion)
