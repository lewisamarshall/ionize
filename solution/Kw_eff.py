def Kw_eff(obj, I=None):
    """Return the effective water dissociation constant.

    Based on the activity corrections to H+ and OH-.
    """
    if not I:
        I = obj.I

    gam_h = obj.H.activity_coefficient(I, 1)
    Kw_eff = obj.Kw/gam_h**2
    return Kw_eff
