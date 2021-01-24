from lensGW.injection.injection_utils import injection_model
import configparser as ConfigParser

class inject(object):
    def __init__(self, config_file):
        self.config_file = config_file
        self.initialize()

    def initialize(self):
        cp = ConfigParser.ConfigParser()
        cp.optionxform = str
        cp.allow_no_value=True
        cp.read(self.config_file)
        self.param = {}  
        #print('----------Param for Noise-----------------\n')
        for (key,val) in cp.items('Param'):
            self.param.update({key: eval(val)})
        #    print(key,':',val)
            
    def do_injection(self,loc_wave,loc_lensed):
        return injection_model(self.param).wave_inject(loc_wave,loc_lensed)
