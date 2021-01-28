import numpy as np
from lensGW.injection.injection import inject
from lensGW.waveform.waveform import gw_signal
import os
import configparser as ConfigParser
import h5py

if __name__ == "__main__" :

    loc_wave = 'lensGW/ini_files/test.ini'
    loc_lensed = 'lensGW/ini_files/lens_param.ini'
    loc_inject = 'lensGW/ini_files/inject.ini'
    
    data_signal = []
    data_lensed_strain = []
    if os.path.isfile(loc_lensed):
        config_file = loc_lensed
        cp = ConfigParser.ConfigParser()
        cp.optionxform = str
        cp.allow_no_value = True
        cp.read(config_file)
        param = {}  
        for (key,val) in cp.items('Param'):
            param.update({key: eval(val)})
        zL = np.linspace(0.2,2,100)
        for i in zL:
            cp.set('Param', 'zL', str(i))
            signal, _ = inject(loc_inject).do_injection(loc_wave,loc_lensed)
            lensed_strain, _, _, _ = gw_signal(loc_wave).lensed_gw(loc_lensed)
            data_signal.append(signal)
            data_lensed_strain.append(lensed_strain)
    else:
        raise RuntimeError('file not found')
        
    hf = h5py.File('data/bd_signal.h5', 'w')
    hf.create_dataset('data', data=data_signal)
    hf.close()
    
    hf = h5py.File('data/bd_strainL.h5', 'w')
    hf.create_dataset('data', data=data_lensed_strain)
    hf.close()

    