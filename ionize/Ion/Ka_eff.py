# def Ka_eff(self, ionic_strength=None, temperature=None):
#     """Return the effective Ka values for the ion.
#
#     Args:
#         I (float): The ambiant ionic strength.
#
#     This function correct the Ka for ionic strength, using the Dubye-Huckel
#     theory to calculate activity coefficients. If no ionic strength is
#     supplied, and the Ion is nested in a Solution, the solution ionic
#     strength will be used. Otherwise, the ionic strength is assumed to be 0.
#     """
#     _, ionic_strength, temperature = \
#         self._resolve_context(None, ionic_strength, temperature)
#
#     # Make the effective Ka vector the same size as the Ka vector.
#     Ka_eff = []
#
#     gam_i = self.activity(ionic_strength)
#     gam_h, = self.activity(1, ionic_strength)
#
#     # For each acidity coefficient, get the effective
#     # coefficient by multiplying by activities.
#     for i, Kp in enumerate(self.Ka(ionic_strength, temperature)):
#         Ka_eff.append(Kp*gam_i[i+1]/gam_i[i]/gam_h)
#
#     return Ka_eff
