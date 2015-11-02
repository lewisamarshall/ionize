"""Constants for use in the Emigrate solver."""
import numpy as np

# Physical Constants
lpm3 = 1000.                    # Liters per meter ** 3
gpkg = 1000.                    # grams per kilogram
faraday = 96485.34              # Faraday's const.[C/mol]
boltzmann = 1.38e-23            # Boltzmann's constant, [J/K]
gas_constant = 8.31             # Universal gas const. [J/mol*K]
permittivity = 8.85e-12         # Permativity of free space. [F/m]
avogadro = 6.02e23              # Avogadro's number, 1/mol
elementary_charge = 1.602e-19   # Charge of a proton, [C]

# Temperature Information
reference_temperature = 25.      # Reference temperature (Celsius)
kelvin_conversion = 273.15

# Correction constants
pitts = 1.5                     # Finite ion radius correction. [(mol/L)**.5]
onsager_fuoss = np.array((0.2929, -0.3536, 0.0884, -0.0442, 0.0276, -0.0193))


def kelvin(temperature):
    """Convert Celsius to Kelvin."""
    return temperature + kelvin_conversion


def celsius(temperature_kelvin):
    """Convert Kelvin to Celsius."""
    return temperature_kelvin - kelvin_conversion

# TODO: Move this information into solvent.
# h_mobility = 362E-9/faraday   # Mobility of Hydroxide   # [m^2/s*V]/F
# oh_mobility = 205E-9/faraday  # Mobility of Hydronium   % [m^2/s*V]/F
# h_diffusivity = h_mobility / 1 * boltzmann * (temperature_K)
# oh_diffusivity = oh_mobility / -1 * boltzmann * (temperature_K)
