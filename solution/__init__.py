classdef solution
    % SOLUTION	Create a  SOLUTION object.
	%
	% 	SOLUTION(IONS, CONCENTRATIONS) is an object representing aqueous solution
	% 		containing the ions in IONS at the concentrations specified in CONCENTRATIONS. 
	%	
	%	When a new solution is initialized, it will immediately calculate the equilibrium
	%		state, including the pH and the ionic strength (I) of the solution. These
	%		values willl be stored as permenant attributes of the object. Other solution 
	%		properties can be calculated by invoking the appropriate method. 
	%	
	% 	See also ion.
    
    properties(Constant = true, GetAccess = 'private')
        F=96485.3415;       % Faraday's const.				[C/mol]
        Rmu=8.31;           % Universal gas const. 			[J/mol*K]
        Temp=298;           % Temperature 					[K]
        Kw=1E-14;           % Water equilibrium constant	[mol^2]
        muH=362E-9;   		% Mobility of Hydronium 		[m^2/s*V]
        muOH=205E-9; 		% Mobility of Hydroxide 		[m^2/s*V]
        Lpm3=1000;          % Liters per meter^3			[]
        visc=1E-3;          % Dynamic viscosity (water) 	[Pa s]
		Adh=0.512; 			% L^1/2 / mol^1/2, approximate for room temperature
		aD=1.5;			 	% mol^-1/2 mol^-3/2, approximation
    end
    
    properties(Constant = false, GetAccess = 'private')
		H=ion('H+', +1, 100, 362E-9);
		OH=ion('OH-', -1, -100, -205E-9);
	end
    properties
        ions;				% Must be a cell of ion objects from the Asp class. 
        concentrations=0; 	% A vector of concentrations in molar.
        pH=7;				% Normal pH units. 
		I=0; 				% Expected in molar. 
    end
    
    methods
        
        function obj=solution(ions, concentrations)
			%Class Constructor

            if(nargin == 2)
				% If the object is not a cell, try to force it to be a cell array. 
				% This may be the case when there is only one ion. 
				if ~iscell(ions)
					ions=num2cell(ions);
				end
				
				% Check that all of the objects in IONS are in the ion class. 
                if isvector(ions) && all(strcmp(cellfun(@class, ions, 'UniformOutput', false), 'ion'))
                    obj.ions=ions;
                else
                    error('You must input a cell vector of ion objects. Use "help ion" for more information.')
                end
				
				% Check that CONCENTRATIONS is a numeric vector of the same length as IONS.
				% Also check that all of the concentrations are positive. 
                if isvector(concentrations) && length(ions)==length(concentrations) && isnumeric(concentrations) && all(concentrations>=0)
					% If the concentration is put in as a cell, change it to a vector. 
					if iscell(concentrations)
						concentrations=cell2mat(concentrations)
					end
					obj.concentrations=concentrations;
                else
                    error('The concentrations vector must be the same size as the ions vector.')
                end
            else % If the solution isn't specified with two arguments, throw an error. 
                error('Solutions must have a cell of ions and a cell or vector of concentrations.')
            end
            
			try
				[obj.pH, obj.I]=obj.find_equilibrium;
			catch
				warning('Could not find equilibrium with ionic strength corrections. Returning uncorrected pH and I. ')
	            obj.pH=obj.calc_pH;
				obj.I=obj.calc_I(obj.pH);
			end

			actual_mobilities=obj.onsager_fuoss;
			for i=1:length(obj.ions)
				obj.ions{i}.actual_mobility=actual_mobilities{i};
			end
			obj.H.actual_mobility=actual_mobilities{end}(1);
			obj.OH.actual_mobility=actual_mobilities{end}(2);
			
        end

		function new_solution=add_ion(obj, new_ions, new_concentrations)
			% ADD_ION initializes a new solution with more ions.
			%	NEW_SOLUTION will contain all of the ions in the current solution
			%	plus a new set of ions from NEW_IONS at a new set of concentrations 
			%	from NEW_CONCENTRATIONS.
			new_solution=solution(cat(2, obj.ions, {new_ions}), cat(2, obj.concentrations, new_concentrations));
		end
				
		function cH=cH(obj, pH, I)
			% Supplies the concentration of H+ in the solution.
			if  ~exist('pH', 'var')
				pH=obj.pH;
			end
			
			if  ~exist('I', 'var')
				I=obj.I;
			end
			
			cH=10^(-pH)/obj.H.activity_coefficient(I, 1);
		end
		
		function cOH=cOH(obj, pH, I)
			% Supplies the concentration of OH- in the solution. 
			if  ~exist('pH', 'var')
				pH=obj.pH;
			end
			if  ~exist('I', 'var')
				I=obj.I;
			end
			
			cOH=obj.Kw_eff(I)/obj.cH(pH);
		end
		
		function H_conductivity=H_conductivity(obj)
			% Calculates teh conductivity of H+.
			% Does not correct the mobility of the ion.
			H_conductivity=obj.cH*obj.H.molar_conductivity(obj.pH, obj.I);
		end
		
		function OH_conductivity=OH_conductivity(obj)
			% Calculates teh conductivity of OH+.
			% Does not correct the mobility of the ion.
			OH_conductivity=obj.cOH*obj.OH.molar_conductivity(obj.pH, obj.I);
		end
        
    end %End methods section
end %End class definition
