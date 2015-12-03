"""Constants used in ionize.

:lpm3:
    1000. liter/m^3.
:gpkg:
    1000. gram/kilogram
:faraday:
    Faraday's constant, 96485. C/mol.
:boltzmann:
    Boltzmann's constant, 1.38e-23 J/K
:gas_constant:
    8.31 J/mol/K
:permittivity:
    Permittivity of free space, 8.854e-12 F/m
:avogadro:
    Avogadro's constant, 6.022e23 / mol.
:elementary_charge:
    Charge of a proton, 1.602e-19 C.
:reference_temperature:
    Room temperature, 25 C.
:pitts: Pitts correction constant for finite ion radius,
    1.5 (mol/L)^.5
"""
import numpy as np

# Physical Constants
lpm3 = 1000.                            # Liters per meter ** 3
gpkg = 1000.                            # grams per kilogram
boltzmann = 1.380649e-23                    # Boltzmann's constant, [J/K]
permittivity = 8.854188e-12                 # Permativity of free space. [F/m]
avogadro = 6.022141e23                      # Avogadro's number, 1/mol
elementary_charge = 1.602177e-19           # Charge of a proton, [C]
gas_constant = boltzmann * avogadro     # Universal gas const. [J/mol*K]
faraday = elementary_charge * avogadro  # Faraday's const.[C/mol]

# Temperature Information
reference_temperature = 25.      # Reference temperature (Celsius)
kelvin_conversion = 273.15

# Environmental Information
atmospheric_CO2 = 0.0004        # Atmospheric CO2, in bar.


# Correction constants
pitts = 1.5                     # Finite ion radius correction. [(mol/L)**.5]
onsager_fuoss = np.array((0.2929, -0.3536, 0.0884, -0.0442, 0.0276, -0.0193))


def kelvin(temperature):
    """Convert Celsius to Kelvin."""
    return temperature + kelvin_conversion


def celsius(temperature_kelvin):
    """Convert Kelvin to Celsius."""
    return temperature_kelvin - kelvin_conversion
