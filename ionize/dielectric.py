def dielectric(obj, T=25):
    """Return the dielectric constant of water at a specified temperature.

    The temperature should be specified in celcius. It is internally converted
    to Kelvin. Correlation is based on the CRC handbook, as implimented in
    STEEP.
    """
    Tk = T+273.15
    d = (2.414e-5)*10**(247.8/(Tk-140))
    return d
