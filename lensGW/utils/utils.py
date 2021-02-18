import numpy as np
import astropy
import sys
import lensGW.constants.constants as const
from lenstronomy.LensModel.lens_model import LensModel

# constants definition
Msun_Kg         = const.M_sun                             # Msun (Kg)
Mpc             = const.Mpc                               # 1 Megaparsec in meters
arcsec2rad      = np.pi/(648000)                          # conversion factor from arcsec to radian
G               = const.G                                 # gravitational constant (m3*Kg-1*s-2)
c               = const.c                                 # speed of light (m/s)

def get_lensed_gws(Fmag, hpt, hct):
        hp_lensed,hc_lensed = np.zeros(len(Fmag),dtype=np.float16),np.zeros(len(Fmag),dtype=np.float16)
        for i in range(len(Fmag)):
            hp_lensed[i], hc_lensed[i] = Fmag[i]*hpt[i], Fmag[i]*hct[i]
        return hp_lensed, hc_lensed
        
def discardOverlaps(inarrX, inarrY, deltas, overlapDist):     
    ToBeRemoved_temp = []
    ToBeRemoved      = []
    
    for i in range(len(inarrX)): 
        item_ra    = inarrX[i]
        item_dec   = inarrY[i]
        
        # consider only the elements after the selected one 
        for j in range(i+1, len(inarrX)): 
            ra  = inarrX[j]
            dec = inarrY[j]

            if np.sqrt(d2([item_ra,item_dec],[ra,dec]))<overlapDist:
                ToBeRemoved_temp.append(j)
            
    # discard index duplicates
    for item in ToBeRemoved_temp:
        if item not in ToBeRemoved:
            ToBeRemoved.append(item)
            
    # finally, delete the overlaps
    inarrX = np.delete(inarrX, ToBeRemoved)
    inarrY = np.delete(inarrY, ToBeRemoved)
    deltas = np.delete(deltas, ToBeRemoved)

    return inarrX, inarrY, deltas
 
def zoom_function(source_pos_x,
                  source_pos_y,
                  grid_width,
                  x_min,
                  y_min,
                  ImgFrS,
                  kwargs_lens,
                  gamma=2,
                  Npixels=30,
                  verbose=False):
    search_window = grid_width*gamma
    min_distance  = search_window/Npixels    
    x_center      = x_min
    y_center      = y_min
    
    res           = ImgFrS.candidate_solutions(sourcePos_x = source_pos_x,
                                               sourcePos_y = source_pos_y,
                                               kwargs_lens         = kwargs_lens,
                                               min_distance        = min_distance,
                                               search_window       = search_window,
                                               x_center            = x_center,
                                               y_center            = y_center,
                                               verbose             = verbose)
    return res
 
 
def magnifications(Img_ra, Img_dec, lens_model_list, kwargs_lens_list, diff=None):
    # instantiate the lens model
    lens_model = LensModel(lens_model_list=lens_model_list)
    
    # magnifications
    mu = lens_model.magnification(Img_ra, Img_dec, kwargs_lens_list, diff=diff)
    
    
    return mu
    

def TimeDelay(Img_ra, Img_dec, source_pos_x, source_pos_y, zL, zS, lens_model_list, kwargs_lens_list, scaled=False, scale_factor=None, cosmo=None): 
    # set a default cosmology if not specified
    if cosmo is None:
        from astropy.cosmology import FlatLambdaCDM              
        cosmo = FlatLambdaCDM(H0=69.7, Om0=0.306, Tcmb0=2.725)
    
    # instantiate the lens model
    lens_model = LensModel(lens_model_list=lens_model_list)
    
    # transform the redshift to angular diameter distance
    DL         = cosmo.angular_diameter_distance(zL)
    DS         = cosmo.angular_diameter_distance(zS)
    DLS        = cosmo.angular_diameter_distance_z1z2(zL, zS)
    D          = DLS/(DL*DS)
    D          = np.float64(D/(Mpc))
    prefactor  = (1+zL)/(D*c)
    shift2     = (Img_ra-source_pos_x)**2+(Img_dec-source_pos_y)**2
    potential  = 0.0

    # compute the time delay, Eq. 6 of Diego et al https://www.aanda.org/articles/aa/pdf/2019/07/aa35490-19.pdf
    potential = lens_model.potential(Img_ra, Img_dec, kwargs_lens_list)    
    td        = prefactor*(0.5*shift2 - potential)
    
    # if scaled, account for it
    if (scaled):
        if scale_factor is None:
            sys.stderr.write('\n\nMust specify a scale factor to use scaled units\n')
            sys.stderr.write('Aborting\n')
            exit(-1)
        else:
            td *= scale_factor**2
            
    # normalize the time delays so that one of the images has delay '0'
    td = td - np.min(td)
        
    return td
    
def NablaTimeDelay(Img_ra, Img_dec, source_pos_x, source_pos_y, zL, zS, lens_model_list, kwargs_lens_list, scaled=False, scale_factor=None, cosmo=None, diff=None): 
    # set a default cosmology if not specified
    if cosmo is None:
        from astropy.cosmology import FlatLambdaCDM              
        cosmo = FlatLambdaCDM(H0=69.7, Om0=0.306, Tcmb0=2.725)
        
    # transform the redshift to angular diameter distance
    DL         = cosmo.angular_diameter_distance(zL)
    DS         = cosmo.angular_diameter_distance(zS)
    DLS        = cosmo.angular_diameter_distance_z1z2(zL, zS)
    D          = DLS/(DL*DS)
    D          = np.float64(D/(Mpc))
    prefactor  = (1+zL)/(D*c)
    shiftX     = Img_ra-source_pos_x
    shiftY     = Img_dec-source_pos_y
    psiX       = 0.0
    psiY       = 0.0
    
    # instantiate the lens model
    lens_model = LensModel(lens_model_list=lens_model_list)
    
    # overall lensing potential derivatives
    psiX, psiY = lens_model.alpha(Img_ra, Img_dec, kwargs_lens_list, diff=diff)
    
    # time delays derivatives, Eq. 6 of Diego et al https://www.aanda.org/articles/aa/pdf/2019/07/aa35490-19.pdf
    tdX = prefactor*(shiftX - psiX)    
    tdY = prefactor*(shiftY - psiY)

    # if scaled, account for it
    if (scaled):
        if scale_factor is None:
            sys.stderr.write('\n\nMust specify a scale factor to use scaled units\n')
            sys.stderr.write('Aborting\n')
            exit(-1)
        else:
            tdX *= scale_factor
            tdY *= scale_factor

    return tdX, tdY
    
def getMinMaxSaddle(Img_ra, Img_dec, lens_model_list, kwargs_lens_list, diff=None):
    # store the results of the test here
    res = []
    
    # instantiate the lens model
    lens_model = LensModel(lens_model_list=lens_model_list)
    
    # second derivative test
    for i in range(len(Img_ra)):
        f_xx, f_xy, f_yx, f_yy = lens_model.hessian(Img_ra[i], Img_dec[i], kwargs_lens_list, diff=diff)
        td_xx = 1-f_xx     # derivatives of the time delay, Eq. (12) in https://iopscience.iop.org/article/10.1086/377430/fulltext/ 
        td_yy = 1-f_yy     
        td_xy = -f_xy
        td_yx = -f_yx
        D = td_xx*td_yy-td_xy*td_yx
        if (D<0):          # it's a saddle point
            res.append(0.5)
        elif (D>0):        
            if (td_xx>0):   # it's a minimum
                res.append(0)
            elif (td_xx<0): # it's a maximum
                res.append(1)
            else:
                sys.stderr.write("\n\nInconclusive test for image # {0}\n".format(i))
                sys.stderr.write("Aborting\n")
                exit(-1)
        else:              # inconclusive test
            if (td_xx>0) or (td_yy>0):
                sys.stderr.write("\n\nInconclusive test for image # {0}. But it is not a maximum\n".format(i))
                sys.stderr.write("Aborting\n")
                exit(-1)
            elif (td_xx<0) or (td_yy<0):
                sys.stderr.write("\n\nInconclusive test for image # {0}. But it is not a minimum\n".format(i))
                sys.stderr.write("Aborting\n")
                exit(-1)
            else:
                sys.stderr.write("\n\nInconclusive test for image # {0}. No possibility can be ruled out.".format(i))
                sys.stderr.write("Aborting\n")
                exit(-1)
                
    return res
 
 
def d2(p0,p1):
    return (p0[0]-p1[0])**2 + (p0[1]-p1[1])**2
    
def param_processing(zL, zS, mL, cosmo=None):
    # set a default flat cosmology if not given
    if cosmo is None: 
        from astropy.cosmology import FlatLambdaCDM
        cosmo = FlatLambdaCDM(H0=69.7, Om0=0.306, Tcmb0=2.725)
    
    # process input parameters    
    DL       = cosmo.angular_diameter_distance(zL)
    DS       = cosmo.angular_diameter_distance(zS)
    DLS      = cosmo.angular_diameter_distance_z1z2(zL, zS)
    D        = DLS/(DL*DS)
    D        = np.float64(D/(Mpc))
    theta_E2 = (4*G*mL*Msun_Kg*D)/c**2                                                                                                                                                                                                                                                                                                                                                                                                                            
    theta_E  = np.sqrt(theta_E2) 
    
    return theta_E