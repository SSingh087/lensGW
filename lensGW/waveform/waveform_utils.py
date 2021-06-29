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
        return self.eval_param(source_ra, source_dec, lens_ra, lens_dec, zS, zL, mL, lens_model_list)
        
    def eval_param(self, source_ra, source_dec, lens_ra, lens_dec, 
                    zS, zL, mL, lens_model_list, optim):
        """
        Finds lensed images for the given set of parameters
        :param source_ra: Right accession of the source of GW (in radians)
        :type source_ra: float
        :param source_dec: Declination of the source of GW (in radians)
        :type source_dec: float
        :param lens_ra: Right accession of the lens (in radians)
        :type lens_ra: array
        :param lens_dec: Declination of the lens (in radians)
        :type lens_dec: array
        :param mL: lens mass
        :type mL: float
        :param zL: lens redshift
        :type zL: float
        :param zS: source redshift
        :type zS: float
        :param lens_model_list: names of the lens profiles to be considered for the lens model
        :type lens_model_list: list of strings
        :param optim: For optimization of search algorithm 
        :type optim: Bool 
        """
    # SECTION BETWEEN --- THIS IS UNDER WORK IN PROGRESS
#--------------------------------------------------------------------------------------------
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

            Img_ra, Img_dec, MacroImg_ra, MacroImg_dec, pixel_width  = microimages(source_pos_x = source_ra,
                                                                                source_pos_y    = source_dec,
                                                                                lens_model_list = lens_model_list,
                                                                                kwargs_lens     = kwargs_lens_list,
                                                                                **solver_kwargs)
#-------------------------------------------------------------------------------------------
        elif len(mL)==1:
            mL, lens_ra, lens_dec = mL[0], lens_ra[0], lens_dec[0]
            thetaE_PM = param_processing(zL, zS, mL)
            kwargs_lens_list = [{'center_x': lens_ra, 'center_y': lens_dec, 'theta_E': thetaE_PM/thetaE_PM}]
            solver_kwargs = {'Scaled'           : True, # indicate that the input is in scaled units 
                             'ScaleFactor'      : thetaE_PM, # and the scale factor  
                             'SearchWindowMacro': 4*thetaE_PM/thetaE_PM,
                             'SearchWindow'     : 4*thetaE_PM/thetaE_PM,
                             'OnlyMacro'        : 'False',
                             'Optimization'     : optim}

            Img_ra, Img_dec, pixel_width = microimages(source_pos_x = source_ra,
                                                                    source_pos_y    = source_dec,
                                                                    lens_model_list = lens_model_list,
                                                                    kwargs_lens     = kwargs_lens_list,
                                                                    **solver_kwargs)
        return Img_ra, Img_dec, kwargs_lens_list, solver_kwargs
