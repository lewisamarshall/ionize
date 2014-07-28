def dielectric(obj, T=25):
    """Return the dielectric constant of water at a specified temperature.

    The temperature should be specified in Celcius. Correlation is based on the
    CRC handbook.
    """
    Tk = T+273.15  # Convert to Kalvin.
    d = 249.21 - 0.79069*Tk + 0.72997e-3*Tk**2
    return d
