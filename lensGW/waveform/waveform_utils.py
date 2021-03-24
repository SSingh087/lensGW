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
            for (key,val) in cp.items('Param'):
                self.param.update({key: eval(val)})

    def param_initialize_and_eval(self):
        y0 = self.param['y0']
        y1 = self.param['y1']
        l0 = self.param['l0']
        l1 = self.param['l1']
        zS = self.param['zS']
        zL = self.param['zL']
        mL = self.param['mL']
        lens_model_list = self.param['lens_model_list']
        return self.eval_param(y0,y1,l0,l1,zS,zL,mL,lens_model_list)
        
    def eval_param(self,y0,y1,l0,l1,zS,zL,mL,lens_model_list,optim='False'):
        if len(mL)>1:
            mtot = sum(mL)
            thetaE  = param_processing(zL, zS, mtot)
            beta0, beta1 = y0*thetaE, y1*thetaE
            thetaE_PM, eta0, eta1 = np.zeros(0), np.zeros(0), np.zeros(0)
            kwargs_lens_list = []
            for i in range(len(mL)):
                thetaE_PM = np.append(thetaE_PM,param_processing(zL, zS, mL[i]))
                eta0 = np.append(eta0,l0[i]*thetaE_PM[i])
                eta1 = np.append(eta1,l1[i]*thetaE_PM[i])
                kwargs_lens_list.append({'center_x': eta0[i],'center_y': eta1[i], 'theta_E': thetaE_PM[i]})
            solver_kwargs = {'SearchWindowMacro': 4*thetaE_PM[0]}
            for i in range(1,len(mL)):
                solver_kwargs.update({'SearchWindow': 4*thetaE_PM[i]})
            solver_kwargs.update({'Optimization': optim})
            Img_ra, Img_dec, MacroImg_ra, MacroImg_dec, pixel_width  = microimages(source_pos_x    = beta0,
                                                                                   source_pos_y    = beta1,
                                                                                   lens_model_list = lens_model_list,
                                                                                   kwargs_lens     = kwargs_lens_list,
                                                                                   **solver_kwargs)
            return Img_ra, Img_dec, beta0, beta1, zL, zS, eta0, eta1, lens_model_list, kwargs_lens_list

        elif len(mL)==1:
            mL,l0,l1 = mL[0],l0[0],l1[0]
            thetaE_PM = param_processing(zL, zS, mL)
            beta0, beta1 = y0*thetaE_PM, y1*thetaE_PM
            eta0, eta1 = l0*thetaE_PM, l1*thetaE_PM
            kwargs_lens_list = [{'center_x': eta0, 'center_y': eta1, 'theta_E': thetaE_PM}]
            solver_kwargs = {'SearchWindowMacro': 4*thetaE_PM,
                             'SearchWindow':      4*thetaE_PM,
                             'Optimization':      optim}
            MacroImg_ra, MacroImg_dec, Macro_pixel_width = microimages(source_pos_x    = beta0,
                                                                        source_pos_y    = beta1,
                                                                        lens_model_list = lens_model_list,
                                                                        kwargs_lens     = kwargs_lens_list,
                                                                        **solver_kwargs)
            return MacroImg_ra, MacroImg_dec, beta0, beta1, zL, zS, eta0, eta1, lens_model_list, kwargs_lens_list
