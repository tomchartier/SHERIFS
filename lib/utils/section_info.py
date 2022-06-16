"""
This modules exctracts some parameters for a given section
"""



def get_nonzero_Mmax(bin_mag,cumulative_rates) :
    """
    gets the maximum magnitude with a non zero rate for this fault.

    bin_mag : list, binning in magnitude for the mfd
    cumulative_rates : cumulative participation rate for the section
    """
    Mmax = bin_mag[0]
    for mag, rate in zip(bin_mag,cumulative_rates):
        if rate != 0.:
            Mmax = mag

    return Mmax
