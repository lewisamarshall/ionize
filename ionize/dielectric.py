def dielectric(obj, T=25):
    """Return the dielectric constant of water at a specified temperature.

    The temperature should be specified in Celcius. Correlation is based on the
    CRC handbook.
    """
    Tk = T+273.15  # Convert to Kalvin.
    d = (2.414e-5)*10**(247.8/(Tk-140))
    return d
