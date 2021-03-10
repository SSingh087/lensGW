import numpy as np
from pycbc.waveform import get_td_waveform, get_fd_waveform
from pycbc.detector import Detector
import configparser as ConfigParser
from lensGW.utils.utils import param_processing
from lensGW.solver.images import microimages

class lens_waveform_model(object):
    
    def __init__(self, config_file):
        if config_file is not None:
            self.config_file = config_file
            cp = ConfigParser.ConfigParser()
            cp.optionxform = str
            cp.allow_no_value=True
            cp.read(self.config_file)
            self.param = {}  
            #print('----------Param for lensed Waveforms-----------------\n')
            for (key,val) in cp.items('Param'):
                self.param.update({key: eval(val)})
                #print(key,':',val)

    def param_initialize_and_eval(self):
        y0 = self.param['y0']
        y1 = self.param['y1']
        l0 = self.param['l0']
        l1 = self.param['l1']
        zS = self.param['zS']
        zL = self.param['zL']
        # masses 
        mL  = self.param['mL']
        lens_model_list = self.param['lens_model_list']
        return self.eval_param(y0,y1,l0,l1,zS,zL,mL,lens_model_list)
        
    def eval_param(self,y0,y1,l0,l1,zS,zL,mL,lens_model_list):
        mL = np.array(mL,dtype=np.float64)
        mtot = mL[0] + mL[1]
        thetaE1 = param_processing(zL, zS, mL[0])
        thetaE2 = param_processing(zL, zS, mL[1])
        thetaE  = param_processing(zL, zS, mtot)

        beta0,beta1 = y0*thetaE,y1*thetaE
        eta10,eta11 = l0*thetaE,l1*thetaE
        eta20,eta21 = -l0*thetaE,l1*thetaE

        kwargs_point_mass_1 = {'center_x': eta10,'center_y': eta11, 'theta_E': thetaE1} 
        kwargs_point_mass_2 = {'center_x': eta20,'center_y': eta21, 'theta_E': thetaE2} 
        kwargs_lens_list    = [kwargs_point_mass_1, kwargs_point_mass_2]  
        solver_kwargs = {'SearchWindowMacro': 4*thetaE1,
                        'SearchWindow'     : 4*thetaE2}
        Img_ra, Img_dec, MacroImg_ra, MacroImg_dec, pixel_width  = microimages(source_pos_x    = beta0,
                                                                                source_pos_y    = beta1,
                                                                                lens_model_list = lens_model_list,
                                                                                kwargs_lens     = kwargs_lens_list,
                                                                                **solver_kwargs)
        #print('Macro Image RA :',MacroImg_ra,'\nMacro Image DEC :',MacroImg_dec,'\npixel width :',pixel_width)
        #print('Solver convergence success !')
        return Img_ra, Img_dec, beta0, beta1, zL, zS, \
                lens_model_list, kwargs_lens_list
