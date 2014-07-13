from scipy.optimize import newton


def find_equilibrium(obj):
    """Return the equilibrium ionic strength and pH.

    It uses the newton's method  root finder from scipy.optimize to find the
    equilibrium pH and ionic strength of a solution, using the ionic-strength
    adjusted activity coefficients. This function is called when the object is
    initialized.
    """
    # Generate an initial ionic strength guess without  activity corrections
    I = obj.calc_I(obj.calc_pH())

    # Iterate to find the true ionic strength.
    I = newton(obj.equil_offset, I)

    # Use this final ionic strength to find the correct pH.
    pH = obj.calc_pH(I)
    return (pH, I)
