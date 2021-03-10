import numpy as np 
from lensGW.utils.utils import TimeDelay, getMinMaxSaddle, magnifications

def amplification_from_data(frequencies, mu, td, n):
    Fmag = np.zeros(len(frequencies))
    for i in range(len(mu)):
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
    
    # will store the results here
    td_list = []   
    mu_list = []
    n_list  = []
    Fmag    = np.zeros(len(frequencies), dtype=np.complex128)
        
    # time delays
    td = np.array(TimeDelay(Img_ra, Img_dec, source_pos_x, source_pos_y, zL, zS, lens_model_list, kwargs_lens_list, scaled=scaled, scale_factor=scale_factor, cosmo=cosmo))

    # magnifications
    mu = np.array(magnifications(Img_ra, Img_dec, lens_model_list, kwargs_lens_list, diff=diff))
    
    # Morse indices
    n  = np.array(getMinMaxSaddle(Img_ra, Img_dec, lens_model_list, kwargs_lens_list, diff=diff))

    # compute the amplification factor
    Fmag = amplification_from_data(frequencies, mu, td, n)

    return Fmag