def Kw_eff(self, I=None):
    """Return the effective water dissociation constant.

    Based on the activity corrections to H+ and OH-.
    """
    if not I:
        I = self.I

    gam_h = self._H.activity_coefficient(I, [1])[0]
    Kw_eff = self._Kw/gam_h**2
    return Kw_eff
