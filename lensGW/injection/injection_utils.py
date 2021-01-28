import pycbc.noise
import pycbc.psd
from pycbc.types import TimeSeries
from lensGW.waveform.waveform import gw_signal
import numpy as np

class injection_model(object):
    def __init__(self,param):
        self.sampling_rate=param['sampling_rate']
        self.flow=param['flow']
        self.delta_f = param['delta_f']
        self.psd = param['psd']
        self.start_time = param['start_time']
        self.delta_t = 1/self.sampling_rate
        self.duration = param['duration']
        self.seed = param['seed']

    def generate_noise(self):
        flen = int((self.sampling_rate/2) / self.delta_f) + 1
        if self.psd == 'aLIGOZeroDetLowPower':
            psd = pycbc.psd.aLIGOZeroDetLowPower(flen, self.delta_f, self.flow)
            tsamples = int(self.duration / self.delta_t)
            self.ts = pycbc.noise.noise_from_psd(tsamples, self.delta_t, psd, self.seed)
            self.ts.start_time = self.start_time
            return self.ts
    
    def wave_inject(self, loc_wave, loc_lensed):
        self.ts = self.generate_noise()
        strain, _,_,_ = gw_signal(loc_wave).unlensed_gw()
        #--------------check if time is in sample times-------------------------#
        
        if strain.sample_times.data[0] in self.ts.sample_times.data and strain.sample_times.data[-1] in self.ts.sample_times.data:
        
        #------------------get lensed waveform-----------------------------------
            lensed_strain, _,_,_ = gw_signal(loc_wave).lensed_gw(loc_lensed)
        
        #---------------------inject here------------------------------------    
            signal, lensed_signal = self.ts.data, self.ts.data
            j=0
            inject_point = np.where(strain.sample_times.data[0] == self.ts.sample_times.data)[0][0]
            for i in range(inject_point,inject_point+len(lensed_strain)): # lensed strain is 1 bin shorter
                signal[i] += strain.data[j]
                lensed_signal[i] += lensed_strain.data[j]
                j+=1
            signal = TimeSeries(signal, delta_t=self.ts.delta_t)
            lensed_signal = TimeSeries(lensed_signal, delta_t=self.ts.delta_t)
            return signal, lensed_signal
        else:
            raise RuntimeError("Check Injection time or length !!")