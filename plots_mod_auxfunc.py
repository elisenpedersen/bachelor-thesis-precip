#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 21 20:00:17 CET 2022


-- Test and plots for mod_auxfunc.F08

@author: Christian Rodehacke, DMI Copenhagen, Denmark

@copyright: Copyright (C) 2016-2022 Christian Rodehacke
"""

# TODO : place all "plotting" subroutines/functions into one common file

# -------------------------------------------------------------------------
def parse_arguments():
    '''
    Parse the command line arguments

    Returns
    -------
    list
        List of strings indicating plot formats: ['png', 'pdf'].

    '''
    parser = argparse.ArgumentParser(
        description='Testing Fortran of CISSEMBEL via c-bindings',
        add_help=False)

    parser.add_argument('--skip_plots', '-s',
                        action='store_false',
                        help='Do not create plots (Default)')

    parser.add_argument('--plot_format', '-f', nargs='?', required=False,
                        default=['png', 'pdf'],
                        choices=['png', 'pdf', 'svg', 'svgz', 'eps', 'ps', 'jpg', 'tif'],
                        help='File format of image output')


    # If --skip_plots or -s is set, set --plot_format=None
    plot_formats = None
    if parser.parse_args().skip_plots:
        plot_formats = parser.parse_args().plot_format

    # # Report used command line arguments
    # for single_arg in vars(parser.parse_args()): #sorted(vars(args)) :
    #     print('  - %-21s =  '% (single_arg)+str(getattr(parser.parse_args(), single_arg)) )

    #If you want to convey all switches outside, use:  return parser.parse_args()
    return plot_formats

# -------------------------------------------------------------------------
def plotting(fig_id, plot_name, suffixes='png', subdirectory='./Figure',
             flag_close=True):
    '''
    Generate plot

    Parameters
    ----------
    fig_id : matplotlib.figure.Figure
        Figure handle.
    plot_name : str
        Prefix of the plot's filename.
    suffixes : str, optional
        Filename suffix defining also the figure file format. The default is 'png'.
    subdirectory : str, optional
        Name of the subdirectory holding the figures. The default is './Figure'.
    flag_close : bool, optional
        Close figure after print, which save memory. The default is True.

    Returns
    -------
    None.

    '''
    print('    Create plot "'+plot_name+'"')

    if not isinstance(suffixes, list):
        suffixes = [suffixes]
    for suffix in suffixes :
        subdir = subdirectory + '/' + suffix.upper()
        if not os.path.exists(subdir) :
            print('      + Create directory '+subdir)
            os.makedirs(subdir)
        else :
            print('      o Reuse directory '+subdir)

        filename = subdir+'/'+plot_name+'.'+suffix
        print('      - Print '+filename)
        fig_id.savefig(filename, orientation='landscape', \
                       transparent='false', facecolor='white') #papertype='a4', \
    if flag_close:
        print('      => Close figure')
        plt.close(fig_id)

# -------------------------------------------------------------------------
def plot_rescale_field(library_name, plot_suffix=None, plotname='rescale_field'):
    '''
    Plot results of the via c-binding call Fortran function subsurface_field4d

    Parameters
    ----------
    library : str or ctypes.CDLL
        Reference to imported shared library containing the C-bindings.
    plot_suffix : str, optional
        Suffix of the figure output file. The default is None.
    plotname : TYPE, optional
        Filename prefix. The default is 'rescale_field'.

    Returns
    -------
    None.

    '''
    print('* Rescale precipitation')
    
    import numpy as np #Elise added this
    
    import cbind_mod_physic as mod_physic #Also added this, but did not remove the error
    
    
    # Inclined plane, shall represent an ice sheet ()
    x = np.minimum(np.linspace(-2000, 4500, 2500), 3500)
    y = np.zeros((1000))
    X, Y = np.meshgrid(x, y)
    del Y

    # Nunatak, dom
    x_nunatak = np.linspace(0, len(x), len(x))
    y_nunatak = np.linspace(0, len(y), len(y))
    X_nunatak, Y_nunatak = np.meshgrid(x_nunatak, y_nunatak)

    # Nunatak x-position left of the center
    x0_nunatak = np.max(x_nunatak)*2./3.+np.min(x_nunatak)*1./3.
    y0_nunatak = np.mean(y_nunatak) # Nunatak y-position at the center

    nunatak = 5000 * \
        np.exp(-np.square(4.0*(X_nunatak-x0_nunatak)/len(x_nunatak))) * \
        np.exp(-np.square(4.0*(Y_nunatak-y0_nunatak)/len(y_nunatak))) - 500.

    # Combined elevation: Inclined plane plus nunatak
    elevation0 = np.maximum(nunatak, np.maximum(X, 0.0))
    elevation_growth = np.maximum(nunatak, np.maximum(X+1000.0, 0.0))
    elevation_decay = np.maximum(nunatak, np.maximum(X-1000.0, 0.0))

    # Build glaciered mask, where
    threshold_elevation2ice = 10
    mask_ice0 = np.where(elevation0 > threshold_elevation2ice, 1, 0)
    mask_ice_growth = np.where(elevation_growth > threshold_elevation2ice, 1, 0)
    mask_ice_decay = np.where(elevation_decay > threshold_elevation2ice, 1, 0)

    # Rainband
    rainband_contour = 1000
    rain_background = np.zeros_like(elevation0) + 5.
    rainband = np.maximum(0., 15.-np.square((elevation0-rainband_contour)*10/len(x)))

    rain = rain_background + rainband

    # Preallocate result fields
    rain_decay = np.zeros_like(rain)
    rain_growth = np.zeros_like(rain)
    rain_decay_rescaled = rain_decay
    rain_growth_rescaled = rain_growth
    mask_decay_rescaled = np.zeros_like(mask_ice0, dtype=int)
    mask_growth_rescaled = np.zeros_like(mask_ice0, dtype=int)

    #
    # Call CISSEMBEL's functions
    #
    use_2dfields = True

    if use_2dfields:
        #
        # Use calls of 2d-fields
        #

        #
        # Height correction of the precipitation
        #
        # Lowered/decayed surface
        rain_decay = cbind_mod_physic.hcprecip_high_desert(rain, elevation0,
                                                     elevation_decay,
                                                     library_name)
        # Uplifted/grown surface
        rain_growth = cbind_mod_physic.hcprecip_high_desert(rain, elevation0,
                                                      elevation_growth,
                                                      library_name)
        #
        # Rescale the precipitation
        #
        # Rescale lowered/decayed surface

        mask_decay_rescaled  = np.maximum(mask_ice0, mask_ice_decay)
        factor = mod_auxfunc.rescale_field(rain, rain_decay, mask_decay_rescaled)
        rain_decay_rescaled = np.where(mask_decay_rescaled,
                                       rain_decay*factor, rain_decay)

        # Rescale uplifted/grown surface
        mask_growth_rescaled = np.maximum(mask_ice0, mask_ice_growth)
        factor = mod_auxfunc.rescale_field(rain, rain_growth, mask_growth_rescaled)

        rain_growth_rescaled = np.where(mask_growth_rescaled,
                                        rain_growth*factor, rain_growth)
    else:
        #
        # Use calls of 1d-fields/arrays
        #

        #
        # Height correction of the precipitation
        #
        for idx in range(rain.shape[0]):
            # Auxillary arrays
            array_rain = rain[idx, :]
            array_elev0 = elevation0[idx, :]
            array_elev_d = elevation_decay[idx, :]
            array_elev_g = elevation_growth[idx, :]

            # Height corrections for precipitation
            # Lowered/decayed surface
            rain_decay[idx, :] = mod_physic.hcprecip_high_desert(array_rain,
                                                                  array_elev0,
                                                                  array_elev_d,
                                                                  library_name)
            # Uplifted/grown surface
            rain_growth[idx, :] = mod_physic.hcprecip_high_desert(array_rain,
                                                                  array_elev0,
                                                                  array_elev_g,
                                                                  library_name)
        del array_elev0, array_elev_d, array_elev_g

        #
        # Rescale the precipitation
        #
        shape_of_elements = rain.shape
        number_of_elements = shape_of_elements[0]*shape_of_elements[1]

        array_rain = rain.reshape(number_of_elements)


        # Rescale lowered/decayed surface
        array_rain_d = rain_decay.reshape(number_of_elements)
        array_mask = np.maximum(mask_ice0.reshape(number_of_elements),
                                mask_ice_decay.reshape(number_of_elements))
        factor = mod_auxfunc.rescale_field(array_rain,
                                           array_rain_d,
                                           array_mask)

        array_rain_decay_rescaled = np.where(array_mask,
                                             array_rain_d*factor,
                                             array_rain_d)

        rain_decay_rescaled = array_rain_decay_rescaled.reshape(shape_of_elements)
        mask_decay_rescaled = array_mask.reshape(shape_of_elements)

        del array_rain_decay_rescaled, array_rain_d

        # Rescale uplifted/grown surface
        array_rain_g = rain_growth.reshape(number_of_elements)
        array_mask = np.maximum(mask_ice0.reshape(number_of_elements),
                                mask_ice_growth.reshape(number_of_elements))

        factor = mod_auxfunc.rescale_field(array_rain,
                                          array_rain_g,
                                          array_mask)

        array_rain_growth_rescaled = np.where(array_mask,
                                              array_rain_g*factor,
                                              array_rain_g)

        rain_growth_rescaled = array_rain_growth_rescaled.reshape(shape_of_elements)
        mask_growth_rescaled = array_mask.reshape(shape_of_elements)

        del array_rain_growth_rescaled, array_rain_g, array_rain, array_mask

    #
    # Figures
    #
    cmap_elevation = 'Greens_r' # 'ocean'
    cmap_rainfall = 'Blues'
    cmap_rainfall_anomaly = 'BrBG'

    color_bandelevation = 'fuchsia'

    linestyle_bandelevation = '-.'

    levels_elevation_contour = [10, 500, 1000, 2000, 3000, 4000]
    levels_elevation_fill = np.arange(0, 5000, 500)
    levels_rainfall = np.arange(4, 24, 2)
    levels_rainfall_anomaly =[-7, -5, -3, -2, -1, 0, 1, 2, 3, 5, 7]

    hatch_no_ice0 = '-'
    hatch_rescaled ='.'

    flag_plot_elevation0 = True
    flag_plot_rainfall0 = True
    flag_plot_elevation_rainfall0 = True
    flag_plot_elevations = True
    flag_plot_rescaled = True
    flag_plot_rescaled_anomaly = True

    def label_nunatak(iaxes, xpos, ypos, print_label=True):
        # Mark Nunatak' peak
        iaxes.plot(xpos, ypos, 'Xk')
        if print_label:
            iaxes.text(xpos, ypos,'  Nunatak', ha='left', va='center')

    def as_si(x, ndp=2):
        # https://stackoverflow.com/questions/31453422/displaying-numbers-with-x-instead-of-e-scientific-notation-in-matplotlib/31453961
        #
        # Call example
        # a=1.92e-7
        # plt.text(0.01, 0.23, r"$a = {0:s}$".format(as_si(a,2)), size=20)

        s = '{x:0.{ndp:d}e}'.format(x=x, ndp=ndp)
        m, e = s.split('e')
        return r'{m:s}\cdot 10^{{{e:d}}}'.format(m=m, e=int(e))


    #
    # Figure 1
    #
    if flag_plot_elevation0:
        fig1 = plt.figure(dpi=300, facecolor='white')
        axes1 = fig1.subplots()
        #axes1.set_aspect('equal')

        # Elevation
        pmesh1 = axes1.contourf(elevation0, cmap=cmap_elevation,
                                levels=levels_elevation_fill)
        cont1 = axes1.contour(elevation0, levels=levels_elevation_contour,
                              colors='dimgray', alpha=0.66, linewidths=2,
                              linestyles='--')

        cont2 = axes1.contour(elevation0, levels=[-999999., rainband_contour],
                              colors='red', alpha=0.66, linewidths=2,
                              linestyles='-')


        # Mask NO ice sheet region
        axes1.contourf(mask_ice0, levels=[0, 0.5, 1],
                       hatches=[hatch_no_ice0, ' '], alpha=0)
        # Mark Nunatak's peak
        label_nunatak(axes1, x0_nunatak, y0_nunatak)

        cbar = fig1.colorbar(pmesh1, ax=axes1)
        cbar.set_label('Elevation (m)')
        if cont1:
            cbar.add_lines(cont1)
        if cont2:
            cbar.add_lines(cont2)

        axes1.xaxis.set_ticks([])
        axes1.yaxis.set_ticks([])
        axes1.text(50, len(y)/2, 'Strips:\nNo ice',
                   color='silver', fontsize='large', va='center')

        axes1.set_title('Reference topography')

        if plot_suffix:
            plotting(fig1, plotname+'_elevation0', plot_suffix)

    #
    # Figure 2
    #
    if flag_plot_rainfall0:
        fig2 = plt.figure(dpi=300, facecolor='white')
        axes1 = fig2.subplots()
        # axes1.set_aspect('equal')

        # Rainfall
        pmesh1 = axes1.contourf(rain, cmap=cmap_rainfall,
                                levels=levels_rainfall, extend='max')
        cont1 = axes1.contour(elevation0, levels=levels_elevation_contour,
                              colors='dimgray', alpha=0.66, linewidths=1,
                              linestyles='--')
        cont2 = axes1.contour(elevation0, levels=[-999999., rainband_contour],
                              colors='red', alpha=0.66, linewidths=2,
                              linestyles='-')

        # Mask NO ice sheet region
        axes1.contourf(mask_ice0, levels=[0, 0.5, 1],
                       hatches=[hatch_no_ice0, ' '], alpha=0)
        # Mark Nunatak's peak
        label_nunatak(axes1, x0_nunatak, y0_nunatak)

        cbar = fig2.colorbar(pmesh1, ax=axes1)
        # cbar.set_label('Rainfall, $\\int p\\; dA=${:4.4g}'.format(rain.sum()))
        cbar.set_label('Rainfall, $\\int p\\; dA={:s}$'.format(as_si(rain.sum())))

        axes1.xaxis.set_ticks([])
        axes1.yaxis.set_ticks([])
        axes1.text(50, len(y)/2, 'Strips:\nNo ice',
                   color='gray', fontsize='large', va='center')

        if plot_suffix:
            plotting(fig2, plotname+'_rainfall0', plot_suffix)


    #
    # Figure 3
    #
    if flag_plot_elevation_rainfall0:
        fig3 = plt.figure(dpi=300, facecolor='white')
        axes1 = fig3.subplots()
        # axes1.set_aspect('equal')

        # Elevation + Rainfall
        # pmesh1 = axes1.pcolormesh(elevation0, cmap=cmap_elevation, vmin=0)
        pmesh1 =axes1.contourf(elevation0, cmap=cmap_elevation,
                               levels=levels_elevation_fill)
        pmesh2 = axes1.pcolormesh(np.where(rain>rain_background, rain, np.NaN),
                                  cmap=cmap_rainfall, vmin=0, alpha=0.22)
        cont1 = axes1.contour(elevation0, levels=levels_elevation_contour,
                              colors='dimgray', alpha=0.66, linewidths=1,
                              linestyles='--')
        cont2 = axes1.contour(elevation0, levels=[-999999., rainband_contour],
                              colors='red', alpha=0.66, linewidths=2,
                              linestyles='-')

        # Mask NO ice sheet region
        axes1.contourf(mask_ice0, levels=[0, 0.5, 1], hatches=[hatch_no_ice0, ' '], alpha=0)
        # Mark Nunatak's peak
        label_nunatak(axes1, x0_nunatak, y0_nunatak)

        cbar = fig3.colorbar(pmesh1, ax=axes1)
        cbar.set_label('Elevation (m)')
        if cont1:
            cbar.add_lines(cont1)
        if cont2:
            cbar.add_lines(cont2)

        axes1.xaxis.set_ticks([])
        axes1.yaxis.set_ticks([])
        axes1.text(50, len(y)/2, 'Strips:\nNo ice',
                   color='silver', fontsize='large', va='center')

        if plot_suffix:
            plotting(fig3, plotname+'elevation0_rain0', plot_suffix)


    #
    # Figure 4
    #
    if flag_plot_elevations:
        fig4 = plt.figure(dpi=300, facecolor='white')
        axes1 = fig4.subplots(3, 1, sharex=True)
        fig4.subplots_adjust(hspace=0.1)

        levels = np.arange(0, 5000, 500)
        # Elevation + Rainfall
        for iax in axes1:
            if iax == axes1[0]:
                pmesh1 = iax.contourf(elevation0, cmap=cmap_elevation,
                                      levels=levels, extend='both')
                cont1 = iax.contour(elevation0,
                                    levels=levels_elevation_contour,
                                    colors='dimgray', alpha=0.66, linewidths=1,
                                    linestyles='--')
                iax.contourf(mask_ice0, levels=[0, 0.5, 1],
                             hatches=[hatch_no_ice0, ' '], alpha=0)
                iax.set_ylabel('Reference')
            elif iax == axes1[1]:
                pmesh2 = iax.contourf(elevation_decay, cmap=cmap_elevation,
                                      levels=levels, extend='both')
                cont1 = iax.contour(elevation_decay,
                                    levels=levels_elevation_contour,
                                    colors='dimgray', alpha=0.66, linewidths=1,
                                    linestyles='--')
                iax.contourf(mask_ice_decay, levels=[0, 0.5, 1],
                             hatches=[hatch_no_ice0, ' '], alpha=0)
                iax.set_ylabel('Shrinking: $\\downarrow$')
            elif iax == axes1[2]:
                pmesh2 = iax.contourf(elevation_growth, cmap=cmap_elevation,
                                      levels=levels, extend='both')
                cont1 = iax.contour(elevation_growth,
                                    levels=levels_elevation_contour,
                                    colors='dimgray', alpha=0.66, linewidths=1,
                                    linestyles='--')
                iax.contourf(mask_ice_growth, levels=[0, 0.5, 1],
                             hatches=[hatch_no_ice0, ' '], alpha=0)
                iax.set_ylabel('Growth: $\\uparrow$')

            cont1.collections[2].set_linestyle(linestyle_bandelevation)
            cont1.collections[2].set_color(color_bandelevation)
            cont1.collections[2].set_linewidth(2)


            cont2 = iax.contour(elevation0,
                                levels=[-999999., rainband_contour],
                                colors='red', alpha=0.66, linewidths=1,
                                linestyles='-')

            # Mark Nunatak's peak
            label_nunatak(iax, x0_nunatak, y0_nunatak)

            iax.xaxis.set_ticks([])
            iax.yaxis.set_ticks([])
            iax.text(50, len(y)/2, 'Strips:\nNo ice',
                     color='silver', fontsize='large', va='center')

        cbar = fig4.colorbar(pmesh1, ax=axes1, location='right')
        cbar.set_label('Elevation (m)')
        if cont1:
            cbar.add_lines(cont1)
        # if cont2:
        #     cbar.add_lines(cont2)

        if plot_suffix:
            plotting(fig4, plotname+'_elevations', plot_suffix)

    #
    # Figure 5
    #
    if flag_plot_rescaled:
        fig5 = plt.figure(dpi=300, facecolor='white')
        axes1 = fig5.subplots(2, 2, sharex=True, sharey= True)
        fig5.subplots_adjust(wspace=0.1)

        for iax in axes1.flatten():
            if iax == axes1[0, 0]:
                pmesh2 = iax.contourf(rain_decay,
                                      cmap=cmap_rainfall,
                                      levels=levels_rainfall, extend='max')
                title = '$\\downarrow\\;\\int p\\, dA={:s}$'.format(as_si(rain_decay.sum()))
                cont1 = iax.contour(elevation_decay,
                                    levels=levels_elevation_contour,
                                    colors='dimgray', alpha=0.66, linewidths=1,
                                    linestyles='--')
            if iax == axes1[0, 1]:
                pmesh2 = iax.contourf(rain_decay_rescaled,
                                      cmap=cmap_rainfall,
                                      levels=levels_rainfall, extend='max')
                title = '$\\downarrow\\Re\\;\\int p\\, dA={:s}$'.format(as_si(rain_decay_rescaled.sum()))
                iax.contourf(mask_decay_rescaled, levels=[0, 0.5, 1],
                              hatches=[' ', hatch_rescaled], alpha=0)
                cont1 = iax.contour(elevation_decay,
                                    levels=levels_elevation_contour,
                                    colors='dimgray', alpha=0.66, linewidths=1,
                                    linestyles='--')


            if iax == axes1[1, 0]:
                pmesh2 = iax.contourf(rain_growth,
                                      cmap=cmap_rainfall,
                                      levels=levels_rainfall, extend='max')
                title = '$\\uparrow\\;\\int p\\, dA={:s}$'.format(as_si(rain_growth.sum()))
                cont1 = iax.contour(elevation_growth,
                                    levels=levels_elevation_contour,
                                    colors='dimgray', alpha=0.66, linewidths=1,
                                    linestyles='--')
            if iax == axes1[1, 1]:
                pmesh2 = iax.contourf(rain_growth_rescaled,
                                      cmap=cmap_rainfall,
                                      levels=levels_rainfall, extend='max')
                title = '$\\uparrow\\Re\\;\\int p\\, dA={:s}$'.format(as_si(rain_growth_rescaled.sum()))
                iax.contourf(mask_growth_rescaled, levels=[0, 0.5, 1],
                              hatches=[' ', hatch_rescaled], alpha=0)
                cont1 = iax.contour(elevation_growth,
                                    levels=levels_elevation_contour,
                                    colors='dimgray', alpha=0.66,
                                    linewidths=2,
                                    linestyles='--')
            if iax == axes1[0, 0]:
                iax.set_ylabel('Shrinking: $\\downarrow$')
            if iax == axes1[1, 0]:
                iax.set_ylabel('Growth $\\uparrow$')
            if iax == axes1[1, 0]:
                iax.set_xlabel('Height correction (HC)')
            if iax == axes1[1, 1]:
                iax.set_xlabel('HC + $\\Re$escaled')

            cont1.collections[2].set_linestyle(linestyle_bandelevation)
            cont1.collections[2].set_color(color_bandelevation)

            cont2 = iax.contour(elevation0,
                                levels=[-999999., rainband_contour],
                                colors='red', alpha=0.66, linewidths=1,
                                linestyles='-')

            iax.set_title(title)

            # Mask NO ice sheet region
            conthatch = iax.contourf(mask_ice0, levels=[0, 0.5, 1],
                         hatches=[hatch_no_ice0, ' '], alpha=0)
            conthatch.collections[1].set_edgecolor('gray')
            # Mark Nunatak's peak
            if iax == axes1[0, 0]:
                label_nunatak(iax, x0_nunatak, y0_nunatak)
            else:
                label_nunatak(iax, x0_nunatak, y0_nunatak, False)

            iax.xaxis.set_ticks([])
            iax.yaxis.set_ticks([])


        cbar = fig5.colorbar(pmesh2, ax=axes1)
        cbar.set_label('Rainfall'+
                       ', Ref: $\\int p\\; dA={:s}$'.format(as_si(rain.sum())))
        # if cont1:
        #     cbar.add_lines(cont1)

        if plot_suffix:
            plotting(fig5, plotname+'_precipitation_rescaled', plot_suffix)


    #
    # Figure 6
    #
    if flag_plot_rescaled_anomaly:
        fig6 = plt.figure(dpi=300, facecolor='white')
        axes1 = fig6.subplots(2, 2, sharex=True, sharey= True)
        fig6.subplots_adjust(wspace=0.1)

        for iax in axes1.flatten():
            if iax == axes1[0, 0]:
                pmesh2 = iax.contourf(rain_decay-rain,
                                      cmap=cmap_rainfall_anomaly,
                                      levels=levels_rainfall_anomaly)
                title = '$\\downarrow\\;\\int p\\, dA={:s}$'.format(as_si(rain_decay.sum()))
                cont1 = iax.contour(elevation_decay,
                                    levels=levels_elevation_contour,
                                    colors='dimgray', alpha=0.66, linewidths=1,
                                    linestyles='--')
            if iax == axes1[0, 1]:
                pmesh2 = iax.contourf(rain_decay_rescaled-rain,
                                      cmap=cmap_rainfall_anomaly,
                                      levels=levels_rainfall_anomaly)
                title = '$\\downarrow\\Re\\;\\int p\\, dA={:s}$'.format(as_si(rain_decay_rescaled.sum()))
                iax.contourf(mask_decay_rescaled, levels=[0, 0.5, 1],
                              hatches=[' ', hatch_rescaled], alpha=0)
                cont1 = iax.contour(elevation_decay,
                                    levels=levels_elevation_contour,
                                    colors='dimgray', alpha=0.66, linewidths=1,
                                    linestyles='--')


            if iax == axes1[1, 0]:
                pmesh2 = iax.contourf(rain_growth-rain,
                                      cmap=cmap_rainfall_anomaly,
                                      levels=levels_rainfall_anomaly)
                title = '$\\uparrow\\;\\int p\\, dA={:s}$'.format(as_si(rain_growth.sum()))
                cont1 = iax.contour(elevation_growth,
                                    levels=levels_elevation_contour,
                                    colors='dimgray', alpha=0.66, linewidths=1,
                                    linestyles='--')
            if iax == axes1[1, 1]:
                pmesh2 = iax.contourf(rain_growth_rescaled-rain,
                                      cmap=cmap_rainfall_anomaly,
                                      levels=levels_rainfall_anomaly)
                title = '$\\uparrow\\Re\\;\\int p\\, dA={:s}$'.format(as_si(rain_growth_rescaled.sum()))
                iax.contourf(mask_growth_rescaled, levels=[0, 0.5, 1],
                              hatches=[' ', hatch_rescaled], alpha=0)
                cont1 = iax.contour(elevation_growth,
                                    levels=levels_elevation_contour,
                                    colors='dimgray', alpha=0.66,
                                    linewidths=2,
                                    linestyles='--')
            if iax == axes1[0, 0]:
                iax.set_ylabel('Shrinking: $\\downarrow$')
            if iax == axes1[1, 0]:
                iax.set_ylabel('Growth $\\uparrow$')
            if iax == axes1[1, 0]:
                iax.set_xlabel('Height correction (HC)')
            if iax == axes1[1, 1]:
                iax.set_xlabel('HC + $\\Re$escaled')

            cont1.collections[2].set_linestyle(linestyle_bandelevation)
            cont1.collections[2].set_color(color_bandelevation)

            cont2 = iax.contour(elevation0,
                                levels=[-999999., rainband_contour],
                                colors='red', alpha=0.66, linewidths=1,
                                linestyles='-')

            iax.set_title(title)

            # Mask NO ice sheet region
            conthatch = iax.contourf(mask_ice0, levels=[0, 0.5, 1],
                         hatches=[hatch_no_ice0, ' '], alpha=0)
            conthatch.collections[1].set_edgecolor('gray')
            # Mark Nunatak's peak
            if iax == axes1[0, 0]:
                label_nunatak(iax, x0_nunatak, y0_nunatak)
            else:
                label_nunatak(iax, x0_nunatak, y0_nunatak, False)

            iax.xaxis.set_ticks([])
            iax.yaxis.set_ticks([])

        cbar = fig6.colorbar(pmesh2, ax=axes1)
        cbar.set_label('Rainfall anomaly'+
                       ', Ref: $\\int p\\; dA={:s}$'.format(as_si(rain.sum())))
        # if cont1:
        #     cbar.add_lines(cont1)

        if plot_suffix:
            plotting(fig6, plotname+'_precipitation_rescaled_ano', plot_suffix)






# -------------------------------------------------------------------------
def plot_values_at_depth_level(library_name, plot_suffix=None,
                               plotname='values_at_depth_level'):
    '''
    Plot results of the via c-binding call Fortran function densprofile

    Parameters
    ----------
    library : str or ctypes.CDLL
        Reference to imported shared library containing the C-bindings.
    plot_suffix : str, optional
        Suffix of the figure output file. The default is None.
    plotname : TYPE, optional
        Filename prefix. The default is 'densprofile'.

    Returns
    -------
    None.

    '''
    print('* Extract values at requested depth')

    depth0 = 0
    depth1 = 50
    depth = np.linspace(depth0, depth1, 11)   # Generic depth profile
    field1 = np.exp(depth/(0.25*max(depth)))  # Generic field1: exponential function
    field2 = np.sin(2*np.pi*depth/max(depth)) # Generic field2: Sin function

    # Higher resolved subsurface depth field which shall reproduce the fields
    # field1 and field2
    depth_subsurface = np.linspace(depth0, depth1, len(depth)*3+2)

    # Preallocte output
    subsurface_values1 = np.zeros_like(depth_subsurface)
    subsurface_values2 = np.zeros_like(depth_subsurface)

    for idepth, sdepth in enumerate(depth_subsurface):
        subsurface_values1[idepth] = \
            mod_auxfunc.values_at_depth_level(field1, depth, sdepth,
                                              library_name)
        subsurface_values2[idepth] = \
            mod_auxfunc.values_at_depth_level(field2, depth, sdepth,
                                             library_name)
    #
    # Figure
    #
    fig = plt.figure(dpi=300, facecolor='white')
    ax_left, ax_right = fig.subplots(1, 2, sharey=True)
    fig.subplots_adjust(wspace=0.025)

    ax_left.plot(field1, depth, '--o', c='cornflowerblue')
    ax_left.plot(subsurface_values1, depth_subsurface, '.', c='darkred')
    ax_left.set_xlabel('Value')
    ax_left.set_title('Exponential', fontsize='small')
    ax_left.legend(['Reference', 'Interpolated'])

    ax_right.plot(field2, depth, '--o', c='cornflowerblue')
    ax_right.plot(subsurface_values2, depth_subsurface, '.', c='darkred')
    ax_right.set_xlabel('Value')
    ax_right.set_title('Sinus', fontsize='small')

    ax_left.set_ylabel('Depth (m)')
    ax_left.invert_yaxis()

    fig.suptitle('Test of interpolation for generic profiles')

    if plot_suffix:
        plotting(fig, plotname, plot_suffix)


# -------------------------------------------------------------------------
#
# Main program
#
if __name__ == '__main__':
    import argparse
    import values_mod_param as physical
    import cbind_mod_physic as mod_physic
    import cbind_mod_auxfunc as mod_auxfunc
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    import os

    print('Load physical constant')
    physical.constant()
    # Values
    GRADLW = physical.constant.gradlw
    LAPSE_RATE = physical.constant.lapse_rate
    PATM = physical.constant.patm_surf
    RHO_ICE = physical.constant.rho_ice
    RHO_SNOW = physical.constant.rho_snow
    RHO_FWATER = physical.constant.rho_fwater
    RHO_SEAWATER = physical.constant.rho_seawater
    CDENSITY = physical.constant.Cdens
    FREEZING_TEMP = physical.constant.Tmelt_fw

    #
    # Load the shared library
    #
    LIBRARY_NAME = '../src/CISSEMBEL_CBind4Tests.so' # INCLUDING path, e.g., ./
    print('Shared library containing C-bindings of Fortran code "'
          +LIBRARY_NAME+'"')

    PLOT_SUFFIXES = parse_arguments() # = None = ['png', 'pdf']

    #
    # Plots for each function
    #

    plot_rescale_field(LIBRARY_NAME, PLOT_SUFFIXES)
    plot_values_at_depth_level(LIBRARY_NAME, PLOT_SUFFIXES)
# -- Last line
