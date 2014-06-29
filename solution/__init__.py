import warnings
from ion import ion


class solution:
    """Represent a solution containing a set of ions.

    Initialize with solution(ions, concentrations), where ions is a list
    containing AspPy ion objects, and concentrations is a list containing
    positive scalars.

    When a new solution is initialized, it will immediately calculate the
    equilibrium state, including the pH and the ionic strength (I) of the
    solution. These values willl be stored as permenant attributes of the
    object. Other solution properties can be calculated by invoking the
    appropriate method.

    See also ion.
    """
    _F = 96485.3415        # Faraday's const.				[C/mol]
    _Rmu = 8.31            # Universal gas const. 	    [J/mol*K]
    _Temp = 298.0          # Temperature 					[K]
    _Kw = 1E-14            # Water equilibrium constant	[mol^2]
    _muH = 362E-9   		 # Mobility of Hydronium 		[m^2/s*V]
    _muOH = 205E-9 		 # Mobility of Hydroxide 		[m^2/s*V]
    _Lpm3 = 1000.0         # Liters per meter^3			[]
    _visc = 1E-3           # Dynamic viscosity (water) 	[Pa s]
    _Adh = 0.512 			 # L^1/2 / mol^1/2, approximate for room temperature
    _aD = 1.5			 	 # mol^-1/2 mol^-3/2, approximation

    _H = ion('H+', +1, 100, 362E-9)
    _OH = ion('OH-', -1, -100, -205E-9)

    ions = []			# Must be a list of ion objects from the Asp class.
    concentrations = 0 	# A list of concentrations in molar.
    pH = 7				# Normal pH units.
    I = 0 				# Expected in molar.

    def __init__(self, ions, concentrations):
        self.ions = ions
        self.concentrations = concentrations
        # Class Constructor

        # if(nargin == 2)
        # 	% If the object is not a cell, try to force it to be a cell array.
        # 	% This may be the case when there is only one ion.
        # 	if ~iscell(ions)
        # 		ions=num2cell(ions);
        # 	end
        #
        # 	% Check that all of the objects in IONS are in the ion class.
        #     if isvector(ions) && all(strcmp(cellfun(@class, ions, 'UniformOutput', false), 'ion'))
        #         obj.ions=ions;
        #     else
        #         error('You must input a cell vector of ion objects. Use "help ion" for more information.')
        #     end
        #
        # 	% Check that CONCENTRATIONS is a numeric vector of the same length as IONS.
        # 	% Also check that all of the concentrations are positive.
        #     if isvector(concentrations) && length(ions)==length(concentrations) && isnumeric(concentrations) && all(concentrations>=0)
        # 		% If the concentration is put in as a cell, change it to a vector.
        # 		if iscell(concentrations)
        # 			concentrations=cell2mat(concentrations)
        # 		end
        # 		obj.concentrations=concentrations;
        #     else
        #         error('The concentrations vector must be the same size as the ions vector.')
        #     end
        # else % If the solution isn't specified with two arguments, throw an error.
        #     error('Solutions must have a cell of ions and a cell or vector of concentrations.')
        # end

        try:
            (obj.pH, obj.I) = obj.find_equilibrium()
        except:
            warnings.warn('''Could not find equilibrium with ionic strength
            corrections. Returning uncorrected pH and I.''')
            obj.pH = obj.calc_pH
            obj.I = obj.calc_I(obj.pH)

        actual_mobilities = obj.onsager_fuoss

        for i in range(len(obj.ions)):
            obj.ions[i].actual_mobility = actual_mobilities[i]

        obj.H.actual_mobility = actual_mobilities[-1][0]
        obj.OH.actual_mobility = actual_mobilities[-1][1]

    def add_ion(obj, new_ions, new_concentrations):
        """add_ion initializes a new solution with more ions.
        NEW_SOLUTION will contain all of the ions in the current solution
        plus a new set of ions from new_ions at a new set of concentrations
        from new_concentrations.
        """
        new_solution = solution((obj.ions + new_ions),
                                [obj.concentrations + new_concentrations])
        return new_solution

    def cH(obj, pH=None, I=None):
        """Return the concentration of H+ in the solution."""
        if not pH:
            pH = obj.pH

        if not I:
            I = obj.I

        cH = 10**(-pH)/obj.H.activity_coefficient(I, 1)[0]
        return cH

    def cOH(obj, pH, I):
        """Supplies the concentration of OH- in the solution."""
        if not pH:
            pH = obj.pH

        if not I:
            I = obj.I

        cOH = obj.Kw_eff(I)/obj.cH(pH)
        return cOH

    def H_conductivity(obj):
        """Calculates teh conductivity of H+.
        Does not correct the mobility of the ion.
        """
        H_conductivity = obj.cH*obj.H.molar_conductivity(obj.pH, obj.I)
        return H_conductivity

    def OH_conductivity(obj):
        """Calculates teh conductivity of OH+.
        Does not correct the mobility of the ion.
        """
        OH_conductivity = obj.cOH*obj.OH.molar_conductivity(obj.pH, obj.I)
        return OH_conductivity
