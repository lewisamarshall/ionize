from scipy.optimize import newton, brentq
import warnings


def find_equilibrium(self):
    """Return the equilibrium ionic strength and pH.

    Uses the newton's method root finder from scipy.optimize to find the
    equilibrium pH and ionic strength of a solution, using the ionic-strength
    adjusted activity coefficients. This function is called when the selfect is
    initialized.
    """
    # Generate an initial ionic strength guess without  activity corrections
    I = self._calc_I(self._calc_pH())

    # Try to bound the true answer.
    b1 = self._equil_offset(0.)
    b2 = self._equil_offset(2.*I)

    # If the answer is in the bound, use brentq. Otherwise, use newton.
    # If brentq doesn't converge, use newton.
    if ((b1 < 0.) ^ (b2 < 0.)):
        I, r = brentq(self._equil_offset, 0, 2*I, full_output=True)
        if not r.converged:
            I = newton(self._equil_offset, I)
            warnings.warn('Couldn\'t use the brentq method. Using newton.')

    else:
        I = newton(self._equil_offset, I)
        warnings.warn('Couldn\'t use the brentq method. Using newton.')

    # Use this final ionic strength to find the correct pH.
    pH = self._calc_pH(I)
    return (pH, I)
