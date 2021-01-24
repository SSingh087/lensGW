import numpy as np
from lensGW.waveform.waveform import gw_signal
from lensGW.injection.injection import inject
from pycbc.types.timeseries import TimeSeries
import h5py

def get_lensed_wave():
    loc_wave = 'lensGW/ini_files/test.ini'
    loc_lensed = 'lensGW/ini_files/lens_param.ini'
    lensed_strain, hp_tilde_lensed, hc_tilde_lensed, freq = gw_signal(loc_wave).lensed_gw(loc_lensed)
    data = [lensed_strain, hp_tilde_lensed, hc_tilde_lensed, freq]
    hf = h5py.File('wave_data.h5', 'w')
    hf.create_dataset('data', data=data)
    hf.close()

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
    print("Creating injections")
    get_injected_wave()
    