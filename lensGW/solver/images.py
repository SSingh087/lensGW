import numpy as np
import sys
from lensGW.solver.solver import lens_eq_solutions
from lenstronomy.LensModel.lens_model import LensModel  
from lensGW.utils.utils import discardOverlaps

def OneDeflector(source_ra,
                 source_dec,
                 lens_model_list,
                 kwargs_lens,
                 **kwargs):

    # see how long it takes to perform the computation
    import time 
    start = time.time() 
    
    # default solverKwargs
    solverKwargs = {'Scaled'           : kwargs['Scaled'],
                    'ScaleFactor'      : kwargs['ScaleFactor'],
                    'SearchWindowMacro': kwargs['SearchWindowMacro'],
                    'PixelsMacro'      : 10**3,
                    'PrecisionLimit'   : 10**(-20),
                    'OverlapDistMacro' : 10**(-15), # prescription: solutions whose distance is less than 10**(-15) rad (\sim 2*1.e-4 micro arcsec) are considered overlaps
                    'NearSource'       : False,
                    'Optimization'     : kwargs['Optimization'],
                    'Verbose'          : False} 
    
    # kwargs that need rescaling if Scaled is True
    ToBeScaled = ['PrecisionLimit', 'OverlapDistMacro']
    
    # update solverKwargs with kwargs
    for key in solverKwargs.keys():
        if key in kwargs:
            value = kwargs[key]
            solverKwargs.update({key: value})
        elif key in ToBeScaled:
            if solverKwargs['Scaled']:
                value = solverKwargs[key]/solverKwargs['ScaleFactor']
                solverKwargs.update({key: value})
        if (key == 'Scaled') and solverKwargs[key]:
            if 'ScaleFactor' not in kwargs:
                sys.stdout.write('\n')
                sys.stdout.write("Scaled units requested, but ScaleFactor not specified: it will be set to {0}\n\n".format(solverKwargs['ScaleFactor']))
    
    if solverKwargs['SearchWindowMacro'] is None:
        sys.stderr.write('\n\nMust specify a search window for the macromodel to perform the analysis\n')
        sys.stderr.write('Aborting\n')
        exit(-1)
     
    # repeat for the optimization-related settings if the optimization mode is active
    if solverKwargs['Optimization']:                           
        optimizationKwargs = {'OptimizationWindowMacro': 2,
                              'OptimizationPixelsMacro': 30,
                              'MinDistMacro': None,
                              'ImprovementMacro': None,
                              'OptimizationPrecisionLimitMacro': 10**(-20)}
                              
        # optimization kwargs that need rescaling if Scaled is True
        ToBeScaledOptimization = ['MinDistMacro', 'OptimizationPrecisionLimitMacro']
        
        # update solverKwargs with optimizationKwargs
        for key in optimizationKwargs.keys():
            if key in kwargs:
                value = kwargs[key]
                solverKwargs.update({key: value})
            elif key in ToBeScaledOptimization:
                if solverKwargs['Scaled']:
                    if optimizationKwargs[key] is None:
                        value = optimizationKwargs[key]
                    else:
                        value = optimizationKwargs[key]/solverKwargs['ScaleFactor']
                    solverKwargs.update({key: value})
                else:
                    value = optimizationKwargs[key]
                    solverKwargs.update({key: value})
            else:
                value = optimizationKwargs[key]
                solverKwargs.update({key: value})
                
    # print them out    
    #sys.stdout.write('\n---- Solver settings ----\n\nThe macromodel analysis will be performed with the following settings:\n\n')
    if solverKwargs['Optimization']:  
        try:
            key = 'OnlyMacro'
            #sys.stdout.write("{0} --> {1}\n".format(key.ljust(len('OptimizationPrecisionLimitMacro')),kwargs[key]))
        except: pass
        for key, value in solverKwargs.items():
            if key != 'ScaleFactor' or (key == 'ScaleFactor' and solverKwargs['Scaled']):
                if key == 'PrecisionLimit'and solverKwargs['Optimization']: pass
                #else:
                    #sys.stdout.write("{0} --> {1}\n".format(key.ljust(len('OptimizationPrecisionLimitMacro')),value))
    else:
        try:
            key = 'OnlyMacro'
            #sys.stdout.write("{0} --> {1}\n".format(key.ljust(len('SearchWindowMacro')),kwargs[key]))
        except: pass
        #for key, value in solverKwargs.items():
            #if key is not 'ScaleFactor' or (key is 'ScaleFactor' and solverKwargs['Scaled']):
                #sys.stdout.write("{0} --> {1}\n".format(key.ljust(len('SearchWindowMacro')),value))
    #sys.stdout.write('\n')

    # parse the remaining relevant arguments 
    search_window_m  = solverKwargs['SearchWindowMacro']    
    min_distance_m   = solverKwargs['SearchWindowMacro']/solverKwargs['PixelsMacro']
        
    # instantiate the lens model
    lens_model  = LensModel(lens_model_list=lens_model_list)
     
    # solve for the images of the macromodel
    #if solverKwargs['Verbose']: sys.stdout.write('\n---- Macromodel analysis ----\n')
    prev_Img_ra, prev_Img_dec, pixel_width = lens_eq_solutions(source_ra  = source_ra,
                                                               source_dec  = source_dec,
                                                               x_center      = source_ra,
                                                               y_center      = source_dec,
                                                               search_window = search_window_m,
                                                               min_distance  = min_distance_m,
                                                               lens_model    = lens_model,
                                                               kwargs_lens   = kwargs_lens,
                                                               macromodel    = True,
                                                               **solverKwargs) 
    # sanity check on the number of solutions
    #sys.stdout.write('\n')
    if len(prev_Img_ra)==0 and solverKwargs['NearSource'] == False:
        sys.stderr.write('\n\nNo images found for the macromodel, try different settings\n')
        sys.stderr.write('Aborting\n')
        exit(-1)
    #elif len(prev_Img_ra) is not 0:   
        #sys.stdout.write('\n\nMACROIMAGES\n\n')
        #sys.stdout.write('ra: {0}\n'.format(prev_Img_ra))
        #sys.stdout.write('dec: {0}\n'.format(prev_Img_dec))
        #sys.stdout.write('\n')
    
    # repeat for the NearSource-related settings if NearSource == True
    if solverKwargs['NearSource']:
        NearSourceKwargs = {'SearchWindowNearSource': None,
                            'PixelsNearSource': 10**3}
                                  
        # update NearSourceKwargs with kwargs
        for key in NearSourceKwargs.keys():
            if key in kwargs:
                value = kwargs[key]
                NearSourceKwargs.update({key: value})
        
        # sanity check on the window near source
        if NearSourceKwargs['SearchWindowNearSource'] is None:
            sys.stderr.write('\n\nMust specify a value for the search window near the source')
            sys.stderr.write('Aborting\n')
            exit(-1)
                
        # update NearSourceKwargs with solverKwargs
        for key in solverKwargs.keys():
            if key != 'SearchWindowMacro' and key != 'PixelsMacro':
                value = solverKwargs[key]
                NearSourceKwargs.update({key: value})
                
        #sys.stdout.write('\n---- Solver settings for the analysis near source ----\n\nThe near source analysis will be performed with the following settings:\n\n')
        #if NearSourceKwargs['Optimization']:  
          #  for key, value in NearSourceKwargs.items():
         #       if key is not 'ScaleFactor' or (key is 'ScaleFactor' and solverKwargs['Scaled']):
                    #sys.stdout.write("{0} --> {1}\n".format(key.ljust(len('OptimizationPrecisionLimitMacro')),value))
        #else:
            #for key, value in NearSourceKwargs.items():
                #if key is not 'ScaleFactor' or (key is 'ScaleFactor' and solverKwargs['Scaled']):
                    #sys.stdout.write("{0} --> {1}\n".format(key.ljust(len('SearchWinodwNearSource')),value))
        sys.stdout.write('\n')     
        
        # parse the remaining relevant arguments
        search_window_ns = NearSourceKwargs['SearchWindowNearSource']
        min_distance_ns  = search_window_ns/NearSourceKwargs['PixelsNearSource']
            
        # check for images very close to the source, which may have been missed when images' displacements are wide
        if solverKwargs['Verbose']: sys.stdout.write('\n---- Near Source analysis ----\n')
        prev_Img_ra_ns, prev_Img_dec_ns, pixel_width_ns = lens_eq_solutions(source_ra  = source_ra,
                                                                            source_dec  = source_dec,
                                                                            x_center      = source_ra,
                                                                            y_center      = source_dec,
                                                                            search_window = search_window_ns,
                                                                            min_distance  = min_distance_ns,
                                                                            lens_model    = lens_model,
                                                                            kwargs_lens   = kwargs_lens,
                                                                            macromodel    = True,
                                                                            **NearSourceKwargs) 
                                                                            
        # append to the previous solutions if there are new ones, then discard the overlaps  
        if len(prev_Img_ra_ns)>0:
        
            #sys.stdout.write('\n')
            #sys.stdout.write('\n\nMACROIMAGES FOUND BY THE NEAR SOURCE CHECK\n\n')
            #sys.stdout.write('ra: {0}\n'.format(prev_Img_ra))
            #sys.stdout.write('dec: {0}\n'.format(prev_Img_dec))
            #sys.stdout.write('\n')
            #sys.stdout.write('\n-----------------------------\n')
            
            prev_Img_ra  = list(prev_Img_ra)
            prev_Img_dec = list(prev_Img_dec)
            
            prev_Img_ra.extend(prev_Img_ra_ns)
            prev_Img_dec.extend(prev_Img_dec_ns)
            
            prev_Img_ra  = np.array(prev_Img_ra)
            prev_Img_dec = np.array(prev_Img_dec)
            
            # not used now, serves as input for the discardOverlaps routine
            dummy_deltas       = np.zeros(len(prev_Img_ra))
            Img_ra, Img_dec, _ = discardOverlaps(prev_Img_ra, prev_Img_dec, dummy_deltas, solverKwargs['OverlapDistMacro'])  
        else:
            #sys.stdout.write('\n')
            #sys.stdout.write('\n\nNO IMAGES FOUND BY THE NEAR SOURCE CHECK\n\n')
            #sys.stdout.write('\n-----------------------------\n\n')
            
            Img_ra  = prev_Img_ra
            Img_dec = prev_Img_dec

    else:
        Img_ra  = prev_Img_ra
        Img_dec = prev_Img_dec
    
    # sanity check on the number of solutions
    if len(Img_ra)==0:
        sys.stderr.write('\n\nNo images found for the macromodel, try different settings\n')
        sys.stderr.write('Aborting\n')
        exit()
    
    #else:
        #if solverKwargs['NearSource']:
            #sys.stdout.write('\n')
            #sys.stdout.write('\n\nTOTAL MACROIMAGES AFTER THE NEAR SOURCE CHECK\n\n')
            #sys.stdout.write('ra: {0}\n'.format(Img_ra))
            #sys.stdout.write('dec: {0}\n'.format(Img_dec))
            #sys.stdout.write('\n')
    
    end = time.time() 
    #if solverKwargs['Verbose']: sys.stdout.write('Elapsed time: {0} seconds\n\n'.format(end-start))

        
    return Img_ra, Img_dec, pixel_width

    
    
def microimages(source_ra,
                source_dec,
                lens_model_list,
                kwargs_lens,
                **kwargs):
    import time 
    
    # default kwargs for the complete model
    solverKwargs = {'Scaled'           : kwargs['Scaled'],  
                    'ScaleFactor'      : kwargs['ScaleFactor'],    
                    'OnlyMacro'        : kwargs['OnlyMacro'], 
                    'MacroIndex'       : [0],                
                    'ImgIndex'         : None,  
                    'SearchWindow'     : kwargs['SearchWindow'],
                    'Pixels'           : 10**3, 
                    'OverlapDist'      : 10**(-15), 
                    'PrecisionLimit'   : 10**(-20), 
                    'Optimization'     : kwargs['Optimization'],
                    'Verbose'          : False} 
    
    # kwargs that need rescaling if Scaled is True
    ToBeScaled = ['OverlapDist', 'PrecisionLimit']
    
    # update the solverkwargs with kwargs
    for key in solverKwargs.keys():
        if key in kwargs:
            value = kwargs[key]
            solverKwargs.update({key: value})
        elif key in ToBeScaled:
            if solverKwargs['Scaled']:
                value = solverKwargs[key]/solverKwargs['ScaleFactor']
                solverKwargs.update({key: value})
                
    # add the extra kwargs
    for key in kwargs.keys():
        if key not in solverKwargs:
            value = kwargs[key]
            solverKwargs.update({key: value})
      
    # set up only the macromodel analysis if all the lenses are defined as macromodel
    if len(lens_model_list) == len(solverKwargs['MacroIndex']):
        solverKwargs.update({'OnlyMacro': True})
        
    # flags for only macromodel/complete analysis
    only_macro  = solverKwargs['OnlyMacro']

    if solverKwargs['SearchWindow'] is None and not only_macro:
            sys.stderr.write('\n\nMust specify a search window for the complete model to perform the analysis\n')
            sys.stderr.write('Aborting\n')
            exit(-1)
        
    #if only_macro:
        #sys.stdout.write('\n---- Will perform only the macromodel analysis ----\n')
        
    # indices of the macromodel components
    macro_index = solverKwargs['MacroIndex']
     
    # solve for the macromodel first
    macromodel_list   = [lens_model_list[index] for index in macro_index]
    macromodel_kwargs = [kwargs_lens[index] for index in macro_index]

    MacroImg_ra, MacroImg_dec, Macro_pixel_width = OneDeflector(source_ra    = source_ra,
                                                                 source_dec    = source_dec,
                                                                 lens_model_list = macromodel_list,
                                                                 kwargs_lens     = macromodel_kwargs,
                                                                 **solverKwargs)   
    
    if only_macro:
        return MacroImg_ra, MacroImg_dec, Macro_pixel_width  
    
    else:
    
        # see how much it takes to perform the computation
        start = time.time()
    
        # settings for the complete model
        model_list    = lens_model_list
        model_kwargs  = kwargs_lens
        search_window = solverKwargs['SearchWindow']
        min_distance  = solverKwargs['SearchWindow']/solverKwargs['Pixels']
        
        # instantiate the complete lens model
        lens_model_complete = LensModel(lens_model_list=model_list)
        
        # index of the macroimage to zoom on
        img_idx       = solverKwargs['ImgIndex']
        
        # update solverKwargs for the complete model analysis 
        if solverKwargs['Optimization']:                           
            optimizationKwargs = {'OptimizationWindow': 2,
                                  'OptimizationPixels': 30,
                                  'MinDist': None,
                                  'Improvement': None,
                                  'OptimizationPrecisionLimit': 10**(-20)}
                                  
            # optimization kwargs that need rescaling if Scaled is True
            ToBeScaledOptimization = ['MinDist', 'OptimizationPrecisionLimit']
            
            # update solverKwargs with optimizationKwargs
            for key in optimizationKwargs.keys():
                if key not in solverKwargs:
                    if key in ToBeScaledOptimization and solverKwargs['Scaled']:
                        if optimizationKwargs[key] is None:
                            value = optimizationKwargs[key]
                        else:
                            value = optimizationKwargs[key]/solverKwargs['ScaleFactor']
                    else:
                        value = optimizationKwargs[key]
                    solverKwargs.update({key: value})

        # print them out
        #sys.stdout.write('\n---- Solver settings ----\n\nThe complete model analysis will be performed with the following settings:\n\n')
        if solverKwargs['Optimization']:  
            for key, value in solverKwargs.items():
                #if key is 'MacroIndex' or key is 'OnlyMacro':
                    #sys.stdout.write("{0} --> {1}\n".format(key.ljust(len('OptimizationPrecisionLimitMacro')),value))
                if 'Macro' in key in key: pass
                elif 'NearSource' in key: pass
                else:
                    if key != 'ScaleFactor' or (key == 'ScaleFactor' and solverKwargs['Scaled']):
                        if key == 'PrecisionLimit' and solverKwargs['Optimization']: pass
                        #else:
                            #sys.stdout.write("{0} --> {1}\n".format(key.ljust(len('OptimizationPrecisionLimitMacro')),value))
        else:
            for key, value in solverKwargs.items():
                #if key is 'MacroIndex' or key is 'OnlyMacro':
                    #sys.stdout.write("{0} --> {1}\n".format(key.ljust(len('PrecisionLimit')),value))
                if 'Macro' in key: pass
                elif 'NearSource' in key: pass
                #else:
                    #if key is not 'ScaleFactor' or (key is 'ScaleFactor' and solverKwargs['Scaled']):
                        #sys.stdout.write("{0} --> {1}\n".format(key.ljust(len('PrecisionLimit')),value))
        #sys.stdout.write('\n')
        
        # pick up the given macroimage, if selected. Iterate over all the macroimages otherwise
        if img_idx  is None:
            img_idx = range(len(MacroImg_ra))
            
        # solve for the images around each macroimage (both for the complete model and the background only)        
        if solverKwargs['Verbose']: sys.stdout.write('\n---- Complete model analysis ----\n\n')
        
        # will store the images before overlap removal
        prev_Img_ra  = [] 
        prev_Img_dec = [] 
        
        # initialize the pixel width to that of the macromodel
        initial_pixel_width = Macro_pixel_width
        
        for idx in img_idx:
            x0, x1       = MacroImg_ra[idx],MacroImg_dec[idx]   
            center_x     = x0
            center_y     = x1

            temp_Img_ra, temp_Img_dec, pixel_width = lens_eq_solutions(source_ra  = source_ra,
                                                                       source_dec  = source_dec,
                                                                       x_center      = center_x,
                                                                       y_center      = center_y,
                                                                       search_window = search_window,
                                                                       min_distance  = min_distance,
                                                                       lens_model    = lens_model_complete,
                                                                       kwargs_lens   = model_kwargs,
                                                                       macromodel    = False,
                                                                       **solverKwargs) 
            
            # store minimum pixel size, in case numerical differentiation is requested
            if pixel_width>initial_pixel_width:
                pixel_width = initial_pixel_width
            initial_pixel_width = pixel_width
        
            # append to the the temporary solutions
            prev_Img_ra.extend(temp_Img_ra)
            prev_Img_dec.extend(temp_Img_dec)

        # discard the overlaps
        dummy_deltas       = np.zeros(len(prev_Img_ra))
        Img_ra, Img_dec, _ = discardOverlaps(prev_Img_ra, prev_Img_dec, dummy_deltas, solverKwargs['OverlapDist'])       

        # sanity check on the number of solutions
        if len(Img_ra)==0:
            sys.stderr.write('\n\nNo images found for the complete model\n')
            sys.stderr.write('Aborting\n')
            exit()
            
        #else:
            #sys.stdout.write('\n')
            #sys.stdout.write('\n\nIMAGES OF THE COMPLETE MODEL\n\n')
            #sys.stdout.write('ra: {0}\n'.format(Img_ra))
            #sys.stdout.write('dec: {0}\n'.format(Img_dec))
            #sys.stdout.write('\n')
            
        end = time.time() 
        if solverKwargs['Verbose']: sys.stdout.write('Elapsed time: {0} seconds\n\n'.format(end-start)) 
        
        return Img_ra, Img_dec, MacroImg_ra, MacroImg_dec, pixel_width      
