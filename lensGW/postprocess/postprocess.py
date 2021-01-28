import numpy as np
import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import os 
import sys
from lenstronomy.LensModel.lens_model import LensModel                                                                                                 
                  
def plot_images(output_folder,
                file_name,
                source_pos_x,
                source_pos_y,
                lens_model_list,
                kwargs_lens_list,
                ImgRA            = None,
                ImgDEC           = None,
                Mu               = None,
                Td               = None,
                MacroImgRA       = None,
                MacroImgDEC      = None,
                MacroMu          = None,
                MacroTd          = None,
                compute_window_x = 10**(-9),
                compute_window_y = 10**(-9),
                center_x         = 0,
                center_y         = 0,
                xmin             = -5*10**(-10),
                ymin             = -4*10**(-10),
                Npixels          = 10**3,
                how_left         = -(0.15*10**(-10)-2*10**(-10)),
                how_right        = 0.15*10**(-10),
                how_up           = 0.2e-11,
                how_down         = 0.7*10**(-10),
                mag_map          = False,
                diff             = None,
                **kwargs):  
    
    ratio = 0.0
    
    # default plotting options
    fig_width_pt  = 3*246.0                                          # get this from LaTeX using \showthe\columnwidth
    inches_per_pt = 1.0/72.27                                        # convert pt to inch
    golden_mean   = ratio if ratio != 0.0 else (np.sqrt(5)-1.0)/2.0  # aesthetic ratio
    fig_width     = fig_width_pt*inches_per_pt                       # width in inches
    fig_height    = fig_width*golden_mean                            # height in inches
    fig_size      = [fig_width,fig_height]

    params = {'axes.labelsize': 22,
              'font.family': 'DejaVu Sans',
              'font.serif': 'Computer Modern Raman',
              'font.size': 20,
              'legend.fontsize': 18,
              'xtick.labelsize': 22,
              'ytick.labelsize': 22,
              'axes.grid' : True,
              'text.usetex': True,
              'savefig.dpi' : 100,
              'lines.markersize' : 14, 
              'axes.formatter.useoffset': False,
              'figure.figsize': fig_size}
    
    # update rcParams with the user defined ones
    try:
        usr_params = kwargs['rcParams']
        for key in usr_params.keys():
            usr_value = usr_params[key]
            params.update({key: usr_value})
    except:
        pass

    mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}'] #for \text command
    mpl.rcParams.update(params)
    
    # parse user defined plotting options
    Macro_color   = kwargs.get('Macro_img_color','dodgerblue')
    Macro_symbol  = kwargs.get('Macro_img_symbol','.')
    img_color      = kwargs.get('img_color','black')
    img_symbol     = kwargs.get('img_symbol','+')
    source_color   = kwargs.get('source_color','red')
    source_symbol  = kwargs.get('source_symbol','.')
    critical_color = kwargs.get('critical_color','black')
    caustics_color = kwargs.get('caustics_color','red')
    title          = kwargs.get('title',None)
    xlabel         = kwargs.get('xlabel',None)
    ylabel         = kwargs.get('ylabel', None)
                                
    # create the output folder if it doesn't already exist
    if output_folder is not None:
        os.system("mkdir -p {0}".format(output_folder))
    
    # instantiate the lens model
    lens_model = LensModel(lens_model_list=lens_model_list)
    
    # compute the caustics and critical curves of the complete model
    from lenstronomy.LensModel.lens_model_extensions import LensModelExtensions

    # maximum resolution => minimum window size
    compute_windowPix   = np.min([compute_window_x, compute_window_y])
    grid_scale          = compute_windowPix/Npixels 
    
    # maximum extension => maximum window size
    compute_window      = np.max([compute_window_x, compute_window_y])  
    
    lensModelExtensions = LensModelExtensions(lens_model)
    ra_crit_list, dec_crit_list,\
    ra_caustic_list, dec_caustic_list = lensModelExtensions.critical_curve_caustics(kwargs_lens_list, compute_window=compute_window, grid_scale=grid_scale)   


    # plot the image positions etc
    plt.plot(source_pos_x, source_pos_y, color=source_color, marker=source_symbol, linestyle='', label=r'$\mathrm{Source}$')
    
    # parse additional datasets, if given and plot them out
    try:
        ras  = kwargs['additional_pt_RAs']
        decs = kwargs['additional_pt_DECs']
        
        try:
            symbols = kwargs['additional_pt_symbols']
        except:
            symbols = ['.' for i in range(len(ras))]
        try:
            colors = kwargs['additional_pt_colors']
        except:
            colors = ['black' for i in range(len(ras))]
        try:
            linestyles = kwargs['additional_pt_linestyles']
        except:
            linestyles = [' ' for i in range(len(ras))]
        try:
            labels = kwargs['additional_pt_labels'] 
        except:
            labels = None
        
        for i in range(len(ras)):
            ra        = ras[i]
            dec       = decs[i]
            symbol    = symbols[i]
            color     = colors[i]
            linestyle = linestyles[i]

            if labels is not None:
                plt.plot(ra, dec, marker=symbol, linestyle=linestyle, color=color, label=labels[i])
            else:
                plt.plot(ra, dec, marker=symbol, linestyle=linestyle, color=color)
    except:
        pass
    
    # Macrourbed images
    if MacroImgRA is not None and MacroImgDEC is not None:
        if len(MacroImgRA)>1:
            plt.plot(MacroImgRA, MacroImgDEC, marker=Macro_symbol, linestyle='', color=Macro_color, label=r'$\mathrm{Macroimages}$')
        else:
            plt.plot(MacroImgRA, MacroImgDEC, marker=Macro_symbol, linestyle='', color=Macro_color, label=r'$\mathrm{Macroimage}$')
        
        # labels  
        if MacroMu is not None:
            if isinstance(MacroMu, int) or isinstance(MacroMu, float):
                MacroMu     = [MacroMu]
                MacroImgRA  = [MacroImgRA]
                MacroImgDEC = [MacroImgDEC]
            for i in range(len(MacroMu)):
                MacroMag =  r'$\mathrm{\mu\,=\,}$'+str(round(np.abs(MacroMu[i]),2))
                if MacroImgRA[i]<center_x:
                    x1 = MacroImgRA[i]-how_left
                else: 
                    x1 = MacroImgRA[i]+how_right
                y1 = MacroImgDEC[i]+how_up
                y2 = MacroImgDEC[i]-how_down
                plt.annotate(MacroMag, xy=(x1, y1))
                
        if MacroTd is not None:
            if isinstance(MacroTd, int) or isinstance(MacroTd, float):
                MacroTd = [MacroTd]
            for i in range(len(MacroTd)):
                MacroTimeDelay = r'$\mathrm{t_d\,=\,}$'+str(round(MacroTd[i],1))+r'$\mathrm{\, d}$'
                if MacroImgRA[i]<center_x:
                    x1 = MacroImgRA[i]-how_left
                else: 
                    x1 = MacroImgRA[i]+how_right
                y1 = MacroImgDEC[i]+how_up
                y2 = MacroImgDEC[i]-how_down
                plt.annotate(MacroTimeDelay, xy=(x1, y2))

    # images of the complete model
    if ImgRA is not None and ImgDEC is not None:
        plt.plot(ImgRA, ImgDEC, marker=img_symbol, mew=1.2, linestyle='', color=img_color, label=r'$\mathrm{Images}$')
    
        # labels 
        if Mu is not None:
            if isinstance(Mu, int) or isinstance(Mu, float):
                Mu = [Mu]
            for i in range(len(Mu)): 
                mag =  r'$\mathrm{\mu\,=\,}$'+str(round(np.abs(Mu[i]),2)) 
                if ImgRA[i]<center_x:
                    x1 = ImgRA[i]-how_left
                else: 
                    x1 = ImgRA[i]+how_right
                y1 = ImgDEC[i]+how_up
                y2 = ImgDEC[i]-how_down
                plt.annotate(mag, xy=(x1, y1))
                
        if Td is not None:
            if isinstance(Td, int) or isinstance(Td, float):
                Td = [Td]
            for i in range(len(Td)):    
                timeDelay = r'$\mathrm{t_d\,=\,}$'+str(round(Td[i]*1.e3,2))+r'$\mathrm{\, ms}$'
                if ImgRA[i]<center_x:
                        x1 = ImgRA[i]-how_left
                else: 
                    x1 = ImgRA[i]+how_right
                y1 = ImgDEC[i]+how_up
                y2 = ImgDEC[i]-how_down      
                plt.annotate(timeDelay, xy=(x1, y2))
        
    # caustics and critical curves
    for i in range(len(ra_crit_list)):
        if i==0:
            plt.plot(ra_crit_list[i], dec_crit_list[i], color=critical_color, linewidth=0.7, label=r'$\mathrm{Critical \,curves}$')
        else:
            plt.plot(ra_crit_list[i], dec_crit_list[i], color=critical_color, linewidth=0.7)
    for i in range(len(ra_caustic_list)):
        if i==0:
            plt.plot(ra_caustic_list[i], dec_caustic_list[i], color=caustics_color, linewidth=0.7, label=r'$\mathrm{Caustics}$')
        else:
            plt.plot(ra_caustic_list[i], dec_caustic_list[i], color=caustics_color, linewidth=0.7)
    
    xmin = center_x-compute_window_x/2. 
    xmax = compute_window_x+xmin
    ymin = center_y-compute_window_y/2.
    ymax = compute_window_y+ymin
    
    # plot out the magnifications
    if mag_map:
        from lensGW.utils.utils import magnifications
        
        num_points_x     = int(compute_window_x/grid_scale)
        num_points_y     = int(compute_window_y/grid_scale)
        num_points       = np.max([num_points_x,num_points_y])
        x,y              = np.linspace(xmin, xmax, num_points), np.linspace(ymin, ymax,num_points)
        x,y              = np.meshgrid(x, y)
        z                = np.zeros((len(x),len(y)))

        for i in range(len(x)):
            for j in range(len(y)):
                mu     = np.abs(magnifications( x[i,j], y[i,j], lens_model_list, kwargs_lens_list, diff=diff))
                z[i,j] = np.abs(np.log10(np.abs(mu)))
        from matplotlib import cm
        plt.contourf(x, y, z, cmap=cm.get_cmap('gray_r'))
        cbar = plt.colorbar()
        cbar.ax.set_ylabel(r'$\mathrm{Magnification \, |\log_{10} \, \mu|}$')
        
    # limits and labels
    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)
    if title is not None:
        plt.title(title)
    if xlabel is not None:
        plt.xlabel(xlabel)
    if ylabel is not None:
        plt.ylabel(ylabel)
    plt.legend(loc='upper left', prop={'size':15})
    plt.tight_layout()
    
    # save or display the plot
    if output_folder is None:
        plt.show()
    else:
        if file_name is None: 
            file_name = 'images.pdf'
        plt.savefig(output_folder+file_name, bbox_inches='tight')
    plt.close()
