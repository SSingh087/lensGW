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
    """
    Computes lensed GW polarizations 
    
    :param F: frequency domain magnification factor :math:`F(f)`, must be computed on the same frequency band as *hptc* and *htc*
    :type F: array
    :param hpt: plus polarization :math:`\\tilde{h}_+(f)`
    :type hpt: array
    :param hct: cross polarization :math:`\\tilde{h}_{\\times}(f)`
    :type hpct: array
    :returns: lensed :math:`\\tilde{h}_+(f), \\tilde{h}_{\\times}(f)` and :math:`s`, computed from :math:`F(f)` and the unlensed quantities
    :rtype: array, array, dict
    """
    hp_lensed = Fmag*hpt.data[:-1]
    hc_lensed = Fmag*hct.data[:-1]

    #frequecy is a bin shorter than the NR waveform 
    #ref https://pycbc.org/pycbc/latest/html/pycbc.waveform.html#pycbc.waveform.utils.frequency_from_polarizations
    return hp_lensed, hc_lensed
        
def discardOverlaps(inarrX, inarrY, deltas, overlapDist):   
    """
    Deletes points in *inarrX* and *inarrY* whose distance is lower than *overlapDist* 
    
    :param inarrX: right ascensions of the points to check (arbitrary units)
    :type inarrX: indexable object
    :param inarrY: declinations of the points to check (arbitrary units)
    :type inarrY: indexable object
    :param deltas: source displacements of the points in *inarrX* and *inarrY* (arbitrary units)
    :type deltas: indexable object
    :param overlapDist: minimum distance for overlap removal
    :type overlapDist: float
    
    :returns: new *inarrX*, *inarrY* and *deltas* which do not contain overlaps
    :rtype: array, array, array
    """
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
    """  
    Zooms to finer grids: computes a new grid from the given specifications and selects pixels in it which could contain solutions of the lens equation
    
    :param source_pos_x: source right ascension (arbitrary units)
    :type source_pos_x: float
    :param source_pos_y: source declination (arbitrary units)
    :type source_pos_y: float
    :param grid_width: grid width before enlargement (arbitrary units)
    :type grid_width: float
    :param x_min: grid center right ascension (arbitrary units)
    :type x_min: float
    :param y_min: grid center right declination (arbitrary units)
    :type y_min: float
    :param ImgFrS: ray-shooting algorithm 
    :type ImgFrS: instance of the ``lenstronomy``'s *LensModel* class
    :param kwargs_lens: keyword arguments of the lens parameters matching each lens profile in the lens model
    :type kwargs_lens: list of dictionaries
    :param gamma: enlargement of the grid width: the grid width will be *grid_width* * *gamma*, *optional*. Default is :math:`2`
    :type gamma: float
    :param Npixels: number of pixels of the new grid, *optional*. Default is :math:`30`
    :type Npixels: int
    :param verbose: prints detailed diagnostic. *Optional*, default is *False*
    :type verbose: bool
      
    :returns: right ascensions, declinations and source displacements of the pixels identified as candidate solutions (arbitrary units) and the pixel width  
    :rtype: array, array, array, float
    """

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
    """     
    Computes image magnifications for a given lens model 
     
    :param Img_ra: images right ascensions (arbitrary units)
    :type Img_ra: array
    :param Img_dec: images declinations (arbitrary units)
    :type Img_dec: array
    :param lens_model_list: names of the lens profiles to be considered for the lens model
    :type lens_model_list: list of strings
    :param kwargs_lens_list: keyword arguments of the lens parameters matching each lens profile in *lens_model_list*
    :type kwargs_lens_list: list of dictionaries
    :param diff: If set, computes the deflection as a finite numerical differential of the lensing
     potential. This differential is only applicable in the single lensing plane where the form of the lensing
     potential is analytically known. *Optional*, default is *None*
    :type diff: None or float
    
    :returns: images magnifications
    :rtype: array
   
    """   
    # instantiate the lens model
    lens_model = LensModel(lens_model_list=lens_model_list)
    
    # magnifications
    mu = lens_model.magnification(Img_ra, Img_dec, kwargs_lens_list, diff=diff)
    
    return mu
    

def TimeDelay(Img_ra, Img_dec, source_pos_x, source_pos_y, zL, zS, lens_model_list, kwargs_lens_list, scaled=False, scale_factor=None, cosmo=None): 
    """    
    Computes image time delays for a given lens model 

    :param Img_ra: images right ascensions (arbitrary units)
    :type Img_ra: array
    :param Img_dec: images declinations (arbitrary units)
    :type Img_dec: array
    :param source_pos_x: source right ascension (arbitrary units)
    :type source_pos_x: float
    :param source_pos_y: source declination (arbitrary units)
    :type source_pos_y: float
    :param zL: lens redshift
    :type zL: float
    :param zS: source redshift
    :type zS: float
    :param lens_model_list: names of the lens profiles to be considered for the lens model
    :type lens_model_list: list of strings
    :param kwargs_lens_list: keyword arguments of the lens parameters matching each lens profile in *lens_model_list*
    :type kwargs_lens_list: list of dictionaries
    :param scaled: specifies if the input is given in arbitrary units, *optional*. If not specified, the input is assumed to be in radians
    :type scaled: bool
    :param scale_factor: scale factor, *optional*. Used to account for the proper conversion factor in the time delays when coordinates are given in arbitrary units, as per :math:`x_{a.u.} = x_{radians}/scale\_factor`. Only considered when *scaled* is *True*
    :type scale_factor: float
    :param cosmo: cosmology used to compute angular diameter distances, *optional*. If not specified, a :math:`\\mathrm{FlatLambdaCDM}` instance with :math:`H_0=69.7, \Omega_0=0.306, T_{cmb0}=2.725` is considered
    :type cosmo: instance of the astropy cosmology class
    
    :returns: time delay in seconds
    :rtype: array
        
    """
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

    """    
    Computes the gradient of the time delay (in seconds) for a given lens model 
    
    :param Img_ra: images right ascensions (arbitrary units)
    :type Img_ra: array
    :param Img_dec: images declinations (arbitrary units)
    :type Img_dec: array
    :param source_pos_x: source right ascension (arbitrary units)
    :type source_pos_x: float
    :param source_pos_y: source declination (arbitrary units)
    :type source_pos_y: float
    :param zL: lens redshift
    :type zL: float
    :param zS: source redshift
    :type zS: float
    :param lens_model_list: names of the lens profiles to be considered for the lens model
    :type lens_model_list: list of strings
    :param kwargs_lens_list: keyword arguments of the lens parameters matching each lens profile in *lens_model_list*
    :type kwargs_lens_list: list of dictionaries
    :param scaled: specifies if the input is given in arbitrary units, *optional*. If not specified, the input is assumed to be in radians
    :type scaled: bool
    :param scale_factor: scale factor, *optional*. Used to account for the proper conversion factor in the time delays when coordinates are given in arbitrary units, as per :math:`x_{a.u.} = x_{radians}/scale\_factor`. Only considered when *scaled* is *True*
    :type scale_factor: float
    :param cosmo: cosmology used to compute angular diameter distances, *optional*. If not specified, a :math:`\\mathrm{FlatLambdaCDM}` instance with :math:`H_0=69.7, \Omega_0=0.306, T_{cmb0}=2.725` is considered
    :type cosmo: instance of the astropy cosmology class
    :param diff: If set, computes the deflection as a finite numerical differential of the lensing
     potential. This differential is only applicable in the single lensing plane where the form of the lensing
     potential is analytically known. *Optional*, default is *None*
    :type diff: None or float
    
    :returns: :math:`\\frac{\\partial t_d}{\\partial x_0},\\frac{\\partial t_d}{\\partial x_1}`
    :rtype: array, array       
    """
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

    """
    Implements the `second derivative test <https://calculus.subwiki.org/wiki/Second_derivative_test_for_a_function_of_two_variables>`_
    to determine if images are minima, maxima or saddle points of the time delay
    
    :param Img_ra: images right ascensions (arbitrary units)
    :type Img_ra: array
    :param Img_dec: images declinations (arbitrary units)
    :type Img_dec: array
    :param lens_model_list: names of the lens profiles to be considered for the lens model
    :type lens_model_list: list of strings
    :param kwargs_lens_list: keyword arguments of the lens parameters matching each lens profile in *lens_model_list*
    :type kwargs_lens_list: list of dictionaries
    :param diff: If set, computes the deflection as a finite numerical differential of the lensing
     potential. This differential is only applicable in the single lensing plane where the form of the lensing
     potential is analytically known. *Optional*, default is *None*
    :type diff: None or float
        
    :returns: list of Morse indices matching the images; :math:`0,1,1/2` if images are minima, maxima or saddle points respectively
    :rtype: list
    
    Raises
    ------
    StandardError
        If the test is inconclusive
    """     
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
    """
    Finds the distance squared between two :math:`2d` points
    
    :param p0: coordinates of the first point
    :type p0: :math:`2d` indexable object
    :param p1: coordinates of the second point
    :type p1: :math:`2d` indexable object
    
    :returns: the distance squared
    :rtype: float
    
    """    
    return (p0[0]-p1[0])**2 + (p0[1]-p1[1])**2
    
def param_processing(zL, zS, mL, cosmo=None):
    """
    Computes the Einstein radius of a given lens
    
    :param zL: lens redshift
    :type zL: float
    :param zS: source redshift
    :type zS: float
    :param mL: lens mass
    :type mL: float
    :param cosmo: cosmology used to compute angular diameter distances, *optional*. If not specified, a :math:`\\mathrm{FlatLambdaCDM}` instance with :math:`H_0=69.7, \Omega_0=0.306, T_{cmb0}=2.725` is considered
    :type cosmo: instance of the astropy cosmology class
    
    :returns: the lens Einstein radius in radians, :math:`\\theta_E = \\sqrt{\\frac{4Gm_L}{c^2}\\frac{D_{LS}}{D_L D_S}}`
    :rtype: float
    """

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