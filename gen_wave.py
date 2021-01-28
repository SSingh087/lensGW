import numpy as np
from lensGW.waveform.waveform import gw_signal
from lensGW.injection.injection import inject
from pycbc.types.timeseries import TimeSeries
import h5py
import pylab

def get_lensed_wave():
    loc_wave = 'lensGW/ini_files/test.ini'
    loc_lensed = 'lensGW/ini_files/lens_param.ini'
    hp_tilde_lensed, hc_tilde_lensed = gw_signal(loc_wave).lensed_gw(loc_lensed)
    lensed_strain = gw_signal(loc_wave).unlensed_strain(loc_lensed, 'H1', 1.13, -1.22, 1.23)
    pylab.plot(lensed_strain.sample_times, lensed_strain)
    pylab.plot(hp_tilde_lensed.sample_times, hp_tilde_lensed)
    pylab.savefig('plot.png')
    
def get_injected_wave():
    loc_wave = 'lensGW/ini_files/test.ini'
    loc_lensed = 'lensGW/ini_files/lens_param.ini'
    loc_inject = 'lensGW/ini_files/inject.ini'
    signal, lensed_signal = inject(loc_inject).do_injection(loc_wave,loc_lensed)
    data = [signal, lensed_signal]
    hf = h5py.File('injection_data.h5', 'w')
    hf.create_dataset('data', data=data)
    hf.close()

if __name__ == "__main__" :
    print("generating waveform")
    get_lensed_wave()