import os
from pycbc import waveform
import configparser as ConfigParser
from lensGW.waveform.waveform_utils import unlens_waveform_model,lens_waveform_model
from lensGW.amplification_factor.amplification_factor import geometricalOpticsMagnification
from lensGW.utils.utils import get_lensed_gws
from pycbc.types.timeseries import TimeSeries

class gw_signal(object):

    def __init__(self, config_file):
        if os.path.isfile(config_file):
            self.config_file = config_file
            self.initialize()
        else:
            raise RuntimeError('file not found')
    
    def initialize(self):
        cp = ConfigParser.ConfigParser()
        cp.optionxform = str
        cp.allow_no_value=True
        cp.read(self.config_file)
        self.param = {}  
        #print('----------Param for Waveforms-----------------\n')
        for (key,val) in cp.items('Param'):
            self.param.update({key: eval(val)})
            #print(key,':',val)
    
    def unlensed_gw(self):
        hp, hc = unlens_waveform_model(self.param).generate()
        print('waveform successfully generated !!')
        return hp, hc
    
    def lensed_gw(self, loc_lensed,
                  diff         = None,
                  scaled       = False,
                  scale_factor = None,
                  cosmo        = None,
                  for_strain   = False):
        if os.path.isfile(loc_lensed):
            hp, hc = self.unlensed_gw()
            freq = waveform.utils.frequency_from_polarizations(hp, hc)
            Img_ra, Img_dec, source_pos_x, source_pos_y,\
            zL, zS, lens_model_list, kwargs_lens_list = lens_waveform_model(loc_lensed).param_initialize_and_eval()
            Fmag = geometricalOpticsMagnification(freq.data,
                                               Img_ra,Img_dec,
                                               source_pos_x,source_pos_y,
                                               zL,zS,
                                               lens_model_list,
                                               kwargs_lens_list,
                                               diff         = diff,
                                               scaled       = scaled,
                                               scale_factor = scale_factor,
                                               cosmo        = cosmo)
            if for_strain: 
                return Fmag
            else :
                #------------return numpy values---------------#
                hp_tilde_lensed, hc_tilde_lensed = get_lensed_gws(Fmag, hp.data, hc.data)
                #------------convert to pycbc.TimeSeries---------------#
                hp_tilde_lensed = TimeSeries(hp_tilde_lensed, delta_t=hp.delta_t)
                hp_tilde_lensed.start_time = hp.start_time
                hc_tilde_lensed = TimeSeries(hc_tilde_lensed, delta_t=hc.delta_t)
                hc_tilde_lensed.start_time = hc.start_time

                return hp_tilde_lensed, hc_tilde_lensed
        else:
            raise RuntimeError('file not found')

    def unlensed_strain(self, loc_lensed, det, ra, dec, polarization):
        from pycbc.detector import Detector
        from numpy import zeros, pad, abs
        hp, hc = hp, hc = self.unlensed_gw()
        strain = Detector(det).project_wave(hp, hc, ra, dec, polarization)
        Fmag = self.lensed_gw(loc_lensed=loc_lensed, for_strain=True)
        Fmag = pad(Fmag, (0, abs(len(Fmag)-len(strain))), 'constant')
        strain_lensed = strain*Fmag #shorter bin here equal to len of Fmag
        return strain_lensed
        
