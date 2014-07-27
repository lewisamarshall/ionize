def viscosity(obj, T=25):
    """Return the viscosity of water at the specified temperature.

    The temperature should be specified in celcius. The function internally
    converts the temperature to Kelvin. Correlation is based on Fox and
    McDonald's Intro to FLuid Mechanics, as implimented in STEEP.
    """
    Tk = T+273.15
    v = (2.414e-5)*10**(247.8/(Tk-140))
    return v
