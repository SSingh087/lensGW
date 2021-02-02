import numpy as np
from pycbc.waveform import get_td_waveform, get_fd_waveform
from pycbc.detector import Detector
import configparser as ConfigParser
from lensGW.utils.utils import param_processing
from lensGW.solver.images import microimages

class unlens_waveform_model(object):

    def __init__(self,param):
        self.domain=param['domain']
        self.approximant=param['approximant']
        self.mass1 = param['mass1']
        self.mass2 = param['mass2']
        self.distance = param['distance']
        self.spin1x = param['spin1x']
        self.spin1y = param['spin1y']
        self.spin1z = param['spin1z']
        self.spin2x = param['spin2x']
        self.spin2y = param['spin2y']
        self.spin2z = param['spin2z']
        self.inclination = param['inclination']
        self.coa_phase = param['coa_phase']
        self.delta_t = param['delta_t']
        self.f_lower = param['f_lower']
        self.end_time = param['end_time']
        self.eccentricity = param['eccentricity']
        self.generate()
    
    def generate(self):
        if self.domain == 'td':
            hp, hc = get_td_waveform(approximant= self.approximant,
                                     mass1= self.mass1,
                                     mass2= self.mass2,
                                     distance=self.distance,
                                     spin1z= self.spin1z,spin1x=self.spin1x,spin1y=self.spin1y,
                                     spin2z= self.spin2z,spin2x=self.spin2x,spin2y=self.spin2y,
                                     inclination= self.inclination,
                                     coa_phase= self.coa_phase,
                                     delta_t= self.delta_t,
                                     f_lower= self.f_lower,
                                     eccentricity = self.eccentricity,
                                    )
        elif self.domain == 'fd':
            hp, hc = get_fd_waveform(approximant= self.approximant,
                                     mass1= self.mass1,
                                     mass2= self.mass2,
                                     distance=self.distance,
                                     spin1z= self.spin1z,spin1x=self.spin1x,spin1y=self.spin1y,
                                     spin2z= self.spin2z,spin2x=self.spin2x,spin2y=self.spin2y,
                                     inclination= self.inclination,
                                     coa_phase= self.coa_phase,
                                     delta_f= 1/self.delta_t,
                                     f_lower= self.f_lower,
                                     eccentricity = self.eccentricity,
                                    )
            hp, hc = hp.to_timeseries(delta_t=self.delta_t), hc.to_timeseries(delta_t=self.delta_t)
        if self.end_time is not None:
            hp.start_time += self.end_time
            hc.start_time += self.end_time

        return hp, hc
            
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
        mtot = mL[0] + mL[1]
        thetaE1 = param_processing(zL, zS, mL[0])
        thetaE2 = param_processing(zL, zS, mL[1])
        thetaE  = param_processing(zL, zS, mtot)

        beta0,beta1 = y0*thetaE,y1*thetaE
        eta10,eta11 = l0*thetaE,l1*thetaE
        eta20,eta21 = -l0*thetaE,l1*thetaE

        # lens model
        
        kwargs_point_mass_1 = {'center_x': eta10,'center_y': eta11, 'theta_E': thetaE1} 
        kwargs_point_mass_2 = {'center_x': eta20,'center_y': eta21, 'theta_E': thetaE2} 
        kwargs_lens_list    = [kwargs_point_mass_1, kwargs_point_mass_2]  

        
        if type(lens_model_list) is list:
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

        else :
            raise RuntimeError("data-type not correct !!")
