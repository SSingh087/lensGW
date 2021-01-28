from lenstronomy.LensModel.Solver.lens_equation_solver import LensEquationSolver
from lensGW.utils.utils import *
import sys
 
def lens_eq_solutions(source_pos_x,
                      source_pos_y,
                      x_center,
                      y_center,
                      search_window,
                      min_distance,
                      lens_model,
                      kwargs_lens,
                      macromodel,
                      **solverKwargs):
    
    # instantiate the ray-shooting class
    ImgFrS = LensEquationSolver(lens_model) 
    
    # initialize the zoom parameters to default values
    gamma   = 2
    Npixels = 30
    
    # extract information for the iterations which is shared between the macromodel and the complete model
    scale_factor = solverKwargs['ScaleFactor']
    verbose      = solverKwargs['Verbose']
    optimization = solverKwargs['Optimization']
        
    # extract information for the iterations which may differ between the macromodel and the complete model
    if macromodel:
        overlap_dist = solverKwargs['OverlapDistMacro']

        if optimization:
            threshold       = solverKwargs['MinDistMacro']  
            improvement     = solverKwargs['ImprovementMacro'] 
            precision_limit = solverKwargs['OptimizationPrecisionLimitMacro']
            try:
                gamma = solverKwargs['OptimizationWindowMacro'] 
            except:
                pass
            try:
                Npixels = solverKwargs['OptimizationPixelsMacro']
            except:
                pass
        else:
            precision_limit = solverKwargs['PrecisionLimit']
    else:
        overlap_dist = solverKwargs['OverlapDist']
        if optimization:
            threshold       = solverKwargs['MinDist']  
            improvement     = solverKwargs['Improvement'] 
            precision_limit = solverKwargs['OptimizationPrecisionLimit']
            try:
                gamma = solverKwargs['OptimizationWindow'] 
            except:
                pass
            try:
                Npixels = solverKwargs['OptimizationPixels']
            except:
                pass
        else:
            precision_limit = solverKwargs['PrecisionLimit']
    
    # find the pixels' centers which minimize locally the displacement (delta_map) w.r.t. the source position
    # returns the pixel width as an input for the iteration    
    x_mins, y_mins, delta_map, pixel_width = ImgFrS.candidate_solutions(sourcePos_x   = source_pos_x,
                                                                        sourcePos_y   = source_pos_y,
                                                                        kwargs_lens   = kwargs_lens,
                                                                        min_distance  = min_distance,
                                                                        search_window = search_window,
                                                                        x_center      = x_center,
                                                                        y_center      = y_center,
                                                                        verbose       = verbose)
                                                                          
    if verbose:
        sys.stdoutz.write('\n')
        sys.stdout.write('Interesting regions of the first grid (no iteration yet):\n\n')
        sys.stdout.write('ra:\n{0}\n'.format(x_mins))
        sys.stdout.write('dec:\n{0}\n'.format(y_mins))
        sys.stdout.write('pixel_width:\n{0}\n'.format(pixel_width))
        sys.stdout.write('\n')

    im_ra           = []    # will store the solutions' ra here
    im_dec          = []    # will store the solutions' dec here
    ctr             = 0     # keeps track of the number of iterations of the subgrid procedure
    check_minsRA    = []    # will contain the ra of the approximate minima around which the zoom is required     
    check_minsDEC   = []    # will contain the dec of the approximate minima around which the zoom is required     
    check_minsDELTA = []    # will contain the ray-shooted distance of the approximate minima around which the zoom is required     
    non_stop        = True  # control variable to stop the iteration if the pixel gets too small
    absmapped_lens  = []    # record the displacements (for DIAGNOSTIC purposes)
    opt_ctr         = 0     # switch on the optimization mode recommendation 

    # check if the pixel size safety threshold has been reached and stop the iteration in that case
    if (pixel_width<10**(-25)/scale_factor):
        #sys.stdout.write('\n\nMinumum pixel size reached. The iteration will be stopped.\n\n')
        non_stop = False 
     
    # iterate on the promising tiles otherwise
    while (non_stop): 
  
        # store solutions of the lens equation which do not require any iteration
        for i in range(len(x_mins)):
            if delta_map[i] < precision_limit: 
                im_ra.extend([x_mins[i]])
                im_dec.extend([y_mins[i]])
                absmapped_lens.extend([delta_map[i]])
            else:
                # minima around which the zoom is required 
                check_minsRA.append(x_mins[i]) 
                check_minsDEC.append(y_mins[i]) 
                check_minsDELTA.append(delta_map[i]) 

        #############################################################
        # optimization procedure, set True for non converging cases #
        #############################################################
        if optimization: 
        
            # discard the pixels which do not satisfy the threshold criterion
            if threshold != None: 

                check_minsRA_temp    = []
                check_minsDEC_temp   = []
                check_minsDELTA_temp = []

                for i in range(len(check_minsRA)): 
                    item_delta = check_minsDELTA[i]

                    if item_delta<threshold:
                        check_minsRA_temp.append(check_minsRA[i])
                        check_minsDEC_temp.append(check_minsDEC[i])
                        check_minsDELTA_temp.append(item_delta)
                           
                check_minsRA    = check_minsRA_temp
                check_minsDEC   = check_minsDEC_temp
                check_minsDELTA = check_minsDELTA_temp

                # update the threshold
                try: threshold = threshold*improvement
                except: pass
       
        ####################################################
        # end of the optimization part                     #
        ####################################################
        
        if (len(check_minsRA)>1000) and (optimization==False) and (opt_ctr==0):
            sys.stdout.write('\n\nMore than 1000 pixels identified as candidate solutions: consider switching to the optimization mode\n\n')
            opt_ctr += 1
            
        if (len(check_minsRA) == 0): 
            # no points on which to iterate
            break
            
        if verbose:
            sys.stdout.write('\n')
            sys.stdout.write('Iteration # {0}\n\n'.format(ctr))
            
            sys.stdout.write('Pixels to iterate over, RAs:')
            sys.stdout.write('\n{0}\n'.format(check_minsRA))
            sys.stdout.write('Pixels to iterate over, DECs:')
            sys.stdout.write('\n{0}\n'.format(check_minsDEC))
            sys.stdout.write('Pixels to iterate over, source displacements:')
            sys.stdout.write('\n{0}\n'.format(check_minsDELTA))
            sys.stdout.write('\n')
                
        # iterate
        ctr+=1 
        
        # don't need the old values anymore. It's time to store the new ones
        x_mins      = [] 
        y_mins      = []
        delta_map   = []

        for k in range(len(check_minsRA)): 
            temp_zoom = zoom_function(source_pos_x = source_pos_x,
                                      source_pos_y = source_pos_y,
                                      grid_width   = pixel_width, 
                                      x_min        = check_minsRA[k],
                                      y_min        = check_minsDEC[k],
                                      ImgFrS       = ImgFrS,
                                      kwargs_lens  = kwargs_lens,
                                      gamma        = gamma,
                                      Npixels      = Npixels,
                                      verbose      = verbose)
                                      
            x_mins.extend(temp_zoom[0])
            y_mins.extend(temp_zoom[1])
            delta_map.extend(temp_zoom[2])
   
        pixel_width = temp_zoom[3] #the same iteration step has the same pixel width
        
        if verbose: sys.stdout.write('\n\nPixel size of the last iteration: {0}\n\n'.format(pixel_width))
        
        if (pixel_width<10**(-25)/scale_factor): 
        
            # stopping condition reached
            #sys.stdout('\n\nMinumum pixel size reached. The iteration will be stopped.\n')
            non_stop = False
            
        else:
        
            # start over again
            check_minsRA    = []                  
            check_minsDEC   = []       
            check_minsDELTA = [] 
    
    # finally, discard the overlaps
    dummy_deltas             = np.zeros(len(im_ra))
    im_ra, im_dec, delta_map = discardOverlaps(im_ra, im_dec, dummy_deltas, overlap_dist)

    return im_ra,im_dec,pixel_width
