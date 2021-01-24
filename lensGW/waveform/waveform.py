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
        print('----------Param for Waveforms-----------------\n')
        for (key,val) in cp.items('Param'):
            self.param.update({key: eval(val)})
            print(key,':',val)
    
    def unlensed_gw(self):
        strain, hp, hc = unlens_waveform_model(self.param).generate()
        print('waveform successfully generated !!')
        freq = waveform.utils.frequency_from_polarizations(hp, hc)
        return strain, hp, hc, freq
    
    def lensed_gw(self, loc_lensed,
                  diff         = None,
                  scaled       = False,
                  scale_factor = None,
                  cosmo        = None):
        if os.path.isfile(loc_lensed):
            strain, hp, hc, freq = self.unlensed_gw()
            Img_ra, Img_dec, source_pos_x, source_pos_y,\
            zL, zS, lens_model_list, kwargs_lens_list, mtot = lens_waveform_model(loc_lensed).param_initialize()
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
            #------------return numpy values---------------#
            hp_tilde_lensed, hc_tilde_lensed, lensed_strain = get_lensed_gws(Fmag, hp.data, hc.data, strain.data)
            #------------convert to pycbc.TimeSeries---------------#
            hp_tilde_lensed = TimeSeries(hp_tilde_lensed, delta_t=hp.delta_t)
            hp_tilde_lensed.start_time = hp.start_time
            hc_tilde_lensed = TimeSeries(hc_tilde_lensed, delta_t=hc.delta_t)
            hc_tilde_lensed.start_time = hc.start_time
            lensed_strain = TimeSeries(lensed_strain, delta_t=strain.delta_t)
            lensed_strain.start_time = strain.start_time
            freq = TimeSeries(Fmag, delta_t=freq.delta_t)
            freq.start_time = freq.start_time
            #waveform.add_custom_waveform('test',waveform_gen,'time', force=True)
            return lensed_strain, hp_tilde_lensed, hc_tilde_lensed, freq
        else:
            raise RuntimeError('file not found')
            
    def waveform_gen(self):
        None