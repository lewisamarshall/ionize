from math import copysign, sqrt
import warnings

from ..constants import pitts


def robinson_stokes_mobility(self, ionic_strength=0., temperature=25.):
    """Return the Robinson-Stokes correction to fully ionized mobility.

    This correction is appropriate if a generic ionic strength is known,
    but the specific ions in solution are unknown.
    """
    # Currently using the ionic strength where Bahga 2010
    # uses twice the ionic strength. This appears to work, and follows the
    # SPRESSO implimentation.
    # Likely typo in paper.
    _, ionic_strength, temperature = \
        self._resolve_context(None, ionic_strength, temperature)

    reference_temperature = self.reference_temperature
    d = self._solvent.dielectric(temperature)
    d_ref = self._solvent.dielectric(reference_temperature)

    A = 0.2297*((reference_temperature+273.15)*d_ref/(temperature+273.15)/d)**(-1.5)
    B = 31.410e-9 * ((reference_temperature+273.15)*d_ref/(temperature+273.15)/d)**(-0.5) *\
        self._solvent.viscosity(reference_temperature)/self._solvent.viscosity(temperature)
    _actual_mobility = []
    for abs_mob, z in zip(self.reference_mobility, self.valence):
        _actual_mobility.append(abs_mob -
                                (A * abs_mob +
                                 B * copysign(1, z)) * sqrt(ionic_strength) /
                                (1. + pitts * sqrt(ionic_strength)))

    return _actual_mobility
