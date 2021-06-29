import numpy as np 
from lensGW.utils.utils import TimeDelay, getMinMaxSaddle, magnifications

def amplification_from_data(frequencies, mu, td, n):
    """
    Computes the geometrical optics amplification :math:`F(f)` over a frequency band, given a set of magnifications, time delays and Morse indices
    
    :param frequencies: frequency band for the computation of the amplification factor
    :type frequencies: array
    :param mu: images magnifications
    :type mu: array
    :param td: images time delays
    :type td: array
    :param n: Morse indices
    :type n: array
    
    :returns: :math:`F(f)`
    :rtype: array
    """
    Fmag = np.zeros(len(frequencies), dtype=np.complex128)
    for i in range(len(mu)):
        #frequecy is a bin shorter than the NR waveform thus so will be Fmag
        #ref https://pycbc.org/pycbc/latest/html/pycbc.waveform.html#pycbc.waveform.utils.frequency_from_polarizations
        Fmag += np.sqrt(np.abs(mu[i]))* np.exp(1j*np.pi*(2.*frequencies*td[i] - n[i]))
    return Fmag

def geometricalOpticsMagnification(frequencies,
                                   Img_ra,
                                   Img_dec,
                                   source_pos_x,
                                   source_pos_y,
                                   zL,
                                   zS,
                                   lens_model_list,
                                   kwargs_lens_list,
                                   diff         = None,
                                   scaled       = False,
                                   scale_factor = None,
                                   cosmo        = None):
    """ Computes the geometrical optics amplification :math:`F(f)` over a frequency band, given a set of images and a lens model
    
    :param frequencies: frequency band for the computation of the amplification factor
    :type frequencies: array
    :param Img_ra: images right ascensions (arbitrary units)
    :type Img_ra: indexable object
    :param Img_dec: images declinations (arbitrary units)
    :type Img_dec: indexable object
    :param source_pos_x: source right ascension (arbitrary units)
    :type source_pos_x: float
    :param source_pos_y: source declination (arbitrary units)
    :type source_pos_y: float
    :param zL: lens redshift
    :type zL: float
    :param zS: source redshift
    :type zS: float
    :param lens_model_list: names of the lens profiles to be considered for the lens model
    :type lens_model_list: list of strings
    :param kwargs_lens_list: keyword arguments of the lens parameters matching each lens profile in *lens_model_list*
    :type kwargs_lens_list: list of dictionaries
    :param diff: step for numerical differentiation, *optional*. Only needed for potentials that require numerical differentiation. If not specified, analytical differentiation is assumed
    :type diff: float
    :param scaled: specifies if the input is given in arbitrary units, *optional*. If not specified, the input is assumed to be in radians
    :type scaled: bool
    :param scale_factor: scale factor, *optional*. Used to account for the proper conversion factor in the time delays when coordinates are given in arbitrary units, as per :math:`x_{a.u.} = x_{radians}/scale\_factor`. Only considered when *scaled* is *True*
    :type scale_factor: float
    :param cosmo: cosmology used to compute angular diameter distances, *optional*. If not specified, a :math:`\\mathrm{FlatLambdaCDM}` instance with :math:`H_0=69.7, \Omega_0=0.306, T_{cmb0}=2.725` is considered
    :type cosmo: instance of the astropy cosmology class
    
    :returns: 
    :math:`F(f)`
    :rtype: array
    """
    
    # will store the results here
    td_list = []   
    mu_list = []
    n_list  = []
    Fmag    = np.zeros(len(frequencies), dtype=np.complex128)
        
    # time delays
    td_list = TimeDelay(Img_ra, Img_dec, source_pos_x, source_pos_y, zL, zS, lens_model_list, kwargs_lens_list, scaled=scaled, scale_factor=scale_factor, cosmo=cosmo)

    # magnifications
    mu_list = magnifications(Img_ra, Img_dec, lens_model_list, kwargs_lens_list, diff=diff)
    
    # Morse indices
    n_list = getMinMaxSaddle(Img_ra, Img_dec, lens_model_list, kwargs_lens_list, diff=diff) 

    td = np.array(td_list)
    mu = np.array(mu_list)
    n  = np.array(n_list)

    # compute the amplification factor
    Fmag = amplification_from_data(frequencies, mu, td, n)

    return Fmag
