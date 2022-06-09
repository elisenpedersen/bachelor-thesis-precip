#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 11:10:19 CEST 2021

C-binding to Fortran functions and subroutines in mod_physic as part of the
Copenhagen Ice Snow Surface Energy and Mass Balance modEL (CISSEMBEL)

@author: Christian Rodehacke, DMI Copenhagen, Denmark
"""
# TODO : Improve python convention, refactor, warning

def densprofile(depth, rho_snow, rho_ice, cdens,
                library='../src/CISSEMBEL_CBind4Tests.so'):
    """
    C-binding to Fortran function densprofile in mod_physic in CISSEMBEL
    Reference density / density profile with depth.

    Parameters
    ----------
    depth : float or numpy.ndarray
        Depth below surface.
    rho_snow : float
        Density of snow.
    rho_ice : float
        Density of ice.
    cdens : float
        Parameter descripting the shape of the density profil.
    library : str or ctypes.CDLL
        Reference to imported shared library containing the C-bindings.

    Returns
    -------
    density : float or numpy.ndarray
        Snow density at given depth(s).

    """
    if isinstance(library, ctypes.CDLL):
        fortlibrary = library
    else:
        fortlibrary = ctypes.cdll.LoadLibrary(library)

    if isinstance(depth, (float, int)):
        inum = 0
        # Single c-value
        c_depth = ctypes.c_double(depth)
    else:
        inum = len(depth)
        # Number conversion, Fortran-order
        c_depth = numpy.array(depth, order="F", dtype=float)
        # c-pointer
        c_depth_ptr = c_depth.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    c_rho_snow = ctypes.c_double(rho_snow)
    c_rho_ice = ctypes.c_double(rho_ice)
    c_cdens = ctypes.c_double(cdens)

    if inum:
        # array
        density = numpy.zeros((inum), order="F", dtype=float)
        density_ptr = density.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        fortlibrary.c_densprofile.argtypes = [ctypes.c_int,
                                              ctypes.POINTER(ctypes.c_double),
                                              ctypes.c_double,
                                              ctypes.c_double,
                                              ctypes.c_double,
                                              ctypes.POINTER(ctypes.c_double)]
        fortlibrary.c_densprofile(ctypes.c_int(inum),
                                  c_depth_ptr,
                                  c_rho_snow,
                                  c_rho_ice,
                                  c_cdens,
                                  density_ptr)
    else:
        # skalar
        fortlibrary.densprofile.restype = ctypes.c_double
        density = fortlibrary.densprofile(ctypes.byref(c_depth),
                                          ctypes.byref(c_rho_ice),
                                          ctypes.byref(c_rho_snow),
                                          ctypes.byref(c_cdens))
    return density


def vaporpress(temp, pres, tmelt0, library='../src/CISSEMBEL_CBind4Tests.so'):
    """
    C-binding to Fortran function vaporpress in mod_physic in CISSEMBEL
    Water vapor pressure.

    Parameters
    ----------
    temp : float or numpy.ndarray
        Temperature.
    pres : float or numpy.ndarray
        Pressure of air.
    tmelt0 : float
        Melting temperature.
    library : str or ctypes.CDLL
        Reference to imported shared library containing the C-bindings.

    Returns
    -------
    water_vapor_pressure : float or numpy.ndarray
        Water vapor pressure at given air pressure level(s).

    """
    if isinstance(library, ctypes.CDLL):
        fortlibrary = library
    else:
        fortlibrary = ctypes.cdll.LoadLibrary(library)

    if isinstance(temp, (float, int)):
        inum = 0
        # Single c-value
        c_temp = ctypes.c_double(temp)
        c_pres = ctypes.c_double(pres)
    else:
        inum = len(temp)
        # Number conversion, Fortran-order
        c_temp = numpy.array(temp, order="F", dtype=float)
        c_pres = numpy.array(pres, order="F", dtype=float)
        # c-pointer
        c_temp_ptr = c_temp.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        c_pres_ptr = c_pres.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    c_tmelt0 = ctypes.c_double(tmelt0)

    if inum:
        # array
        water_vapor_pressure = numpy.zeros((inum), order="F", dtype=float)
        water_vapor_pressure_ptr = \
            water_vapor_pressure.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        fortlibrary.c_vaporpress.argtypes = [ctypes.c_int,
                                             ctypes.POINTER(ctypes.c_double),
                                             ctypes.POINTER(ctypes.c_double),
                                             ctypes.c_double,
                                             ctypes.POINTER(ctypes.c_double)]
        fortlibrary.c_vaporpress(ctypes.c_int(inum),
                                 c_temp_ptr,
                                 c_pres_ptr,
                                 c_tmelt0,
                                 water_vapor_pressure_ptr)
    else:
        # skalar
        fortlibrary.vaporpress.restype = ctypes.c_double
        water_vapor_pressure = fortlibrary.vaporpress(ctypes.byref(c_temp),
                                                      ctypes.byref(c_pres),
                                                      ctypes.byref(c_tmelt0))

    return water_vapor_pressure


def rhoair(pres, library='../src/CISSEMBEL_CBind4Tests.so'):
    """
    C-binding to Fortran function rhoair in mod_physic in CISSEMBEL
    Density of air.

    Parameters
    ----------
    pres : float or numpy.ndarray
        Pressure of air.
    library : str or ctypes.CDLL
        Reference to imported shared library containing the C-bindings.

    Returns
    -------
    air_density : float or numpy.ndarray
        Density of air at given pressure level(s).

    """
    if isinstance(library, ctypes.CDLL):
        fortlibrary = library
    else:
        fortlibrary = ctypes.cdll.LoadLibrary(library)

    if isinstance(pres, (float, int)):
        inum = 0
        # Single c-value
        c_pres = ctypes.c_double(pres)
    else:
        inum = len(pres)
        # Number conversion, Fortran-order
        c_pres = numpy.array(pres, order="F", dtype=float)
        # c-pointer
        c_pres_ptr = c_pres.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

    if inum:
        # array
        air_density = numpy.zeros((inum), order="F", dtype=float)
        air_density_ptr = air_density.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        fortlibrary.c_rhoair.argtypes = [ctypes.c_int,
                                         ctypes.POINTER(ctypes.c_double),
                                         ctypes.POINTER(ctypes.c_double)]
        fortlibrary.c_rhoair(ctypes.c_int(inum), c_pres_ptr, air_density_ptr)
    else:
        # skalar
        fortlibrary.rhoair.restype = ctypes.c_double
        air_density = fortlibrary.rhoair(ctypes.byref(c_pres))
    return air_density

def temp_excess2melt(dtemp, density, thick, cp_layer,
                     library='../src/CISSEMBEL_CBind4Tests.so'):
    """
    C-binding to Fortran function texcess2melt in mod_physic in CISSEMBEL
    Temperature excess leads to melting.

    Parameters
    ----------
    dtemp : float or numpy.ndarray
        Temperature difference.
    density : float or numpy.ndarray
        Layer's density.
    thick : float or numpy.ndarray
        Layer's thickness.
    cp_layer : float or numpy.ndarray
        Layer's heat capacity.
    library : str or ctypes.CDLL
        Reference to imported shared library containing the C-bindings.

    Returns
    -------
    melting : float or numpy.ndarray
        Melting.

    """
    if isinstance(library, ctypes.CDLL):
        fortlibrary = library
    else:
        fortlibrary = ctypes.cdll.LoadLibrary(library)

    if isinstance(dtemp, float) and False: # TODO todo : add skalar-function to mod_physic
        inum = 0
        # Single c-value
        c_dtemp = ctypes.c_double(dtemp)
        c_density = ctypes.c_double(density)
        c_thick = ctypes.c_double(thick)
        c_cp_layer = ctypes.c_double(cp_layer)
    else:
        inum = len(dtemp)
        # Number conversion, Fortran-order
        c_dtemp = numpy.array(dtemp, order="F", dtype=float)
        c_density = numpy.array(density, order="F", dtype=float)
        c_thick = numpy.array(thick, order="F", dtype=float)
        c_cp_layer = numpy.array(cp_layer, order="F", dtype=float)
        # c-pointer
        c_dtemp_ptr = c_dtemp.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        c_density_ptr = c_density.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        c_thick_ptr = c_thick.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        c_cp_layer_ptr = c_cp_layer.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

    if inum:
        # array
        melting = numpy.zeros((inum), order="F", dtype=float)
        melting_ptr = \
            melting.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        fortlibrary.c_texcess2melt.argtypes = [ctypes.c_int,
                                               ctypes.POINTER(ctypes.c_double),
                                               ctypes.POINTER(ctypes.c_double),
                                               ctypes.POINTER(ctypes.c_double),
                                               ctypes.POINTER(ctypes.c_double),
                                               ctypes.POINTER(ctypes.c_double)]
        fortlibrary.c_texcess2melt(ctypes.c_int(inum),
                                   c_dtemp_ptr,
                                   c_density_ptr,
                                   c_thick_ptr,
                                   c_cp_layer_ptr,
                                   melting_ptr)
    else:
        # skalar
        fortlibrary.temp_excess2melt.restype = ctypes.c_double
        melting = fortlibrary.texcess2melt(ctypes.byref(c_dtemp),
                                           ctypes.byref(c_density),
                                                      ctypes.byref(c_thick),
                                                      ctypes.byref(c_cp_layer))
    return melting


def sensheatflux(rho, wind, tair, tsuf, library='../src/CISSEMBEL_CBind4Tests.so'):
    """
    C-binding to Fortran function sensheatflux in mod_physic in CISSEMBEL
    Sensible heat flux.

    Parameters
    ----------
    rho : float or numpy.ndarray
        Density of air.
    wind : float or numpy.ndarray
        Absolute wind velocity.
    tair : float or numpy.ndarray
        Air temperature.
    tsuf : float or numpy.ndarray
        Surface temperature.
    library : str or ctypes.CDLL
        Reference to imported shared library containing the C-bindings.

    Returns
    -------
    sensible_heat_flux : float or numpy.ndarray
        Sensible heat flux.

    """
    if isinstance(library, ctypes.CDLL):
        fortlibrary = library
    else:
        fortlibrary = ctypes.cdll.LoadLibrary(library)

    if isinstance(rho, float) and False: # TODO todo : add skalar-function to mod_physic
        inum = 0
        # Single c-value
        c_rho = ctypes.c_double(rho)
        c_wind = ctypes.c_double(wind)
        c_tair = ctypes.c_double(tair)
        c_tsuf = ctypes.c_double(tsuf)
    else:
        inum = len(rho)
        # Number conversion, Fortran-order
        c_rho = numpy.array(rho, order="F", dtype=float)
        c_wind = numpy.array(wind, order="F", dtype=float)
        c_tair = numpy.array(tair, order="F", dtype=float)
        c_tsuf = numpy.array(tsuf, order="F", dtype=float)
        # c-pointer
        c_rho_ptr = c_rho.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        c_wind_ptr = c_wind.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        c_tair_ptr = c_tair.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        c_tsuf_ptr = c_tsuf.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

    if inum:
        # array
        sensible_heat_flux = numpy.zeros((inum), order="F", dtype=float)
        sensible_heat_flux_ptr = \
            sensible_heat_flux.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        fortlibrary.c_sensheatflux.argtypes = [ctypes.c_int,
                                               ctypes.POINTER(ctypes.c_double),
                                               ctypes.POINTER(ctypes.c_double),
                                               ctypes.POINTER(ctypes.c_double),
                                               ctypes.POINTER(ctypes.c_double),
                                               ctypes.POINTER(ctypes.c_double)]
        fortlibrary.c_sensheatflux(ctypes.c_int(inum),
                                   c_rho_ptr,
                                   c_wind_ptr,
                                   c_tair_ptr,
                                   c_tsuf_ptr,
                                   sensible_heat_flux_ptr)
    else:
        # skalar
        fortlibrary.sensheatflux.restype = ctypes.c_double
        sensible_heat_flux = fortlibrary.sensheatflux(ctypes.byref(c_rho),
                                                      ctypes.byref(c_wind),
                                                      ctypes.byref(c_tair),
                                                      ctypes.byref(c_tsuf))
    return sensible_heat_flux


def latheatflux(wind, vap_air, vap_suf, library='../src/CISSEMBEL_CBind4Tests.so'):
    """
    C-binding to Fortran function latheatflux in mod_physic in CISSEMBEL
    Latent heat flux.

    Parameters
    ----------
    wind : float or numpy.ndarray
        Absolute wind velocity.
    vap_air : float or numpy.ndarray
        Water vapor pressure in the air.
    vap_suf : float or numpy.ndarray
        Water vapor pressure at the surface.
    library : str or ctypes.CDLL
        Reference to imported shared library containing the C-bindings.

    Returns
    -------
    latent_heat_flux : float or numpy.ndarray
        Latent heat flux.

    """
    if isinstance(library, ctypes.CDLL):
        fortlibrary = library
    else:
        fortlibrary = ctypes.cdll.LoadLibrary(library)


    if isinstance(wind, float) and False: # TODO todo : add skalar-function to mod_physic
        inum = 0
        # Single c-value
        c_wind = ctypes.c_double(wind)
        c_vap_air = ctypes.c_double(vap_air)
        c_vap_suf = ctypes.c_double(vap_suf)
    else:
        inum = len(wind)
        # Number conversion, Fortran-order
        c_wind = numpy.array(wind, order="F", dtype=float)
        c_vap_air = numpy.array(vap_air, order="F", dtype=float)
        c_vap_suf = numpy.array(vap_suf, order="F", dtype=float)
        # c-pointer
        c_wind_ptr = c_wind.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        c_vap_air_ptr = c_vap_air.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        c_vap_suf_ptr = c_vap_suf.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

    if inum:
        # array
        latent_heat_flux = numpy.zeros((inum), order="F", dtype=float)
        latent_heat_flux_ptr = \
            latent_heat_flux.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        fortlibrary.c_latheatflux.argtypes = [ctypes.c_int,
                                              ctypes.POINTER(ctypes.c_double),
                                              ctypes.POINTER(ctypes.c_double),
                                              ctypes.POINTER(ctypes.c_double),
                                              ctypes.POINTER(ctypes.c_double)]
        fortlibrary.c_latheatflux(ctypes.c_int(inum),
                                  c_wind_ptr,
                                  c_vap_air_ptr,
                                  c_vap_suf_ptr,
                                  latent_heat_flux_ptr)
    else:
        # skalar
        fortlibrary.latheatflux.restype = ctypes.c_double
        latent_heat_flux = fortlibrary.latheatflux(ctypes.byref(c_wind),
                                                   ctypes.byref(c_vap_air),
                                                   ctypes.byref(c_vap_suf))
    return latent_heat_flux


def sublimation(lhflux, stemp, mtemp, library='../src/CISSEMBEL_CBind4Tests.so'):
    """
    C-binding to Fortran function sublimation in mod_physic in CISSEMBEL
    Sublimation.

    Parameters
    ----------
    lhflux : float or numpy.ndarray
        Latent heat flux.
    stemp : float or numpy.ndarray
        Surface temperature.
    mtemp : float or numpy.ndarray
        Melting temperature.
    library : str or ctypes.CDLL
        Reference to imported shared library containing the C-bindings.

    Returns
    -------
    sublimation_flux : float or numpy.ndarray
        Sublimation (inkl depostion aka negative sublimation).

    """
    if isinstance(library, ctypes.CDLL):
        fortlibrary = library
    else:
        fortlibrary = ctypes.cdll.LoadLibrary(library)

    if isinstance(lhflux, float) and False: # TODO todo : add skalar-function to mod_physic
        inum = 0
        # Single c-value
        c_lhflux = ctypes.c_double(lhflux)
        c_stemp = ctypes.c_double(stemp)
        c_mtemp = ctypes.c_double(mtemp)
    else:
        inum = len(lhflux)
        # Number conversion, Fortran-order
        c_lhflux = numpy.array(lhflux, order="F", dtype=float)
        c_stemp = numpy.array(stemp, order="F", dtype=float)
        c_mtemp = numpy.array(mtemp, order="F", dtype=float)
        # c-pointer
        c_lhflux_ptr = c_lhflux.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        c_stemp_ptr = c_stemp.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        c_mtemp_ptr = c_mtemp.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

    if inum:
        # array
        sublimation_flux = numpy.zeros((inum), order="F", dtype=float)
        sublimation_flux_ptr = \
            sublimation_flux.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        fortlibrary.c_sublimation_flux.argtypes = [ctypes.c_int,
                                              ctypes.POINTER(ctypes.c_double),
                                              ctypes.POINTER(ctypes.c_double),
                                              ctypes.POINTER(ctypes.c_double),
                                              ctypes.POINTER(ctypes.c_double)]
        fortlibrary.c_sublimation_flux(ctypes.c_int(inum),
                                       c_lhflux_ptr,
                                       c_stemp_ptr,
                                       c_mtemp_ptr,
                                       sublimation_flux_ptr)
    else:
        # skalar
        fortlibrary.sublimation_flux.restype = ctypes.c_double
        sublimation_flux = fortlibrary.sublimation_flux(ctypes.byref(c_lhflux),
                                              ctypes.byref(c_stemp),
                                              ctypes.byref(c_mtemp))
    return sublimation_flux


def evaporation(lhflux, stemp, mtemp, library='../src/CISSEMBEL_CBind4Tests.so'):
    """
    C-binding to Fortran function evaporation in mod_physic in CISSEMBEL
    Evaporation.

    Parameters
    ----------
    lhflux : float or numpy.ndarray
        Latent heat flux.
    stemp : float or numpy.ndarray
        Surface temperature.
    mtemp : float or numpy.ndarray
        Melting temperature.
    library : str or ctypes.CDLL
        Reference to imported shared library containing the C-bindings.

    Returns
    -------
    evaporation_flux : float or numpy.ndarray
        Evaporation (inkl. condensation aka negative evaporation).

    """
    if isinstance(library, ctypes.CDLL):
        fortlibrary = library
    else:
        fortlibrary = ctypes.cdll.LoadLibrary(library)

    if isinstance(lhflux, float) and False: # TODO todo : add skalar-function to mod_physic
        inum = 0
        # Single c-value
        c_lhflux = ctypes.c_double(lhflux)
        c_stemp = ctypes.c_double(stemp)
        c_mtemp = ctypes.c_double(mtemp)
    else:
        inum = len(lhflux)
        # Number conversion, Fortran-order
        c_lhflux = numpy.array(lhflux, order="F", dtype=float)
        c_stemp = numpy.array(stemp, order="F", dtype=float)
        c_mtemp = numpy.array(mtemp, order="F", dtype=float)
        # c-pointer
        c_lhflux_ptr = c_lhflux.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        c_stemp_ptr = c_stemp.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        c_mtemp_ptr = c_mtemp.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

    if inum:
        # array
        evaporation_flux = numpy.zeros((inum), order="F", dtype=float)
        evaporation_flux_ptr = \
            evaporation_flux.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        fortlibrary.c_evaporation_flux.argtypes = [ctypes.c_int,
                                              ctypes.POINTER(ctypes.c_double),
                                              ctypes.POINTER(ctypes.c_double),
                                              ctypes.POINTER(ctypes.c_double),
                                              ctypes.POINTER(ctypes.c_double)]
        fortlibrary.c_evaporation_flux(ctypes.c_int(inum),
                                       c_lhflux_ptr,
                                       c_stemp_ptr,
                                       c_mtemp_ptr,
                                       evaporation_flux_ptr)
    else:
        # skalar
        fortlibrary.evaporation_flux.restype = ctypes.c_double
        evaporation_flux = fortlibrary.evaporation_flux(ctypes.byref(c_lhflux),
                                              ctypes.byref(c_stemp),
                                              ctypes.byref(c_mtemp))
    return evaporation_flux


def tpre2snowf_cut(tprec, tair, tmelt0, library='../src/CISSEMBEL_CBind4Tests.so'):
    """
    C-binding to Fortran function tpre2snowf_cut in mod_physic in CISSEMBEL
    Snowfall rate of total precipitation based on air temperature threshold.

    Parameters
    ----------
    tprec : float or numpy.ndarray
        Total precipitation.
    tair : float or numpy.ndarray
        Air temperature.
    tmelt0 : float
        Melting temperature separating snowfall and rainfall.
    library : str or ctypes.CDLL
        Reference to imported shared library containing the C-bindings.

    Returns
    -------
    tprecip2snowf_cut : float or numpy.ndarray
        Snowfall rate from total precipitation.

    """
    if isinstance(library, ctypes.CDLL):
        fortlibrary = library
    else:
        fortlibrary = ctypes.cdll.LoadLibrary(library)

    if isinstance(tprec, (float, int)) and False:
        inum = 0
        # Single c-value
        c_tprec = ctypes.c_double(tprec)
        c_tair = ctypes.c_double(tair)
    else:
        inum = len(tprec)
        # Number conversion, Fortran-order
        c_tprec = numpy.array(tprec, order="F", dtype=float)
        c_tair = numpy.array(tair, order="F", dtype=float)
        # c-pointer
        c_tprec_ptr = c_tprec.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        c_tair_ptr = c_tair.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    c_tmelt0 = ctypes.c_double(tmelt0)

    if inum:
        # array
        tprecip2snowf_cut = numpy.zeros((inum), order="F", dtype=float)
        tprecip2snowf_cut_ptr = \
            tprecip2snowf_cut.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        fortlibrary.c_tpre2snowf_cut.argtypes = [ctypes.c_int,
                                                 ctypes.POINTER(ctypes.c_double),
                                                 ctypes.POINTER(ctypes.c_double),
                                                 ctypes.c_double,
                                                 ctypes.POINTER(ctypes.c_double)]
        fortlibrary.c_tpre2snowf_cut(ctypes.c_int(inum),
                                     c_tprec_ptr,
                                     c_tair_ptr,
                                     c_tmelt0,
                                     tprecip2snowf_cut_ptr)
    else:
        # skalar
        fortlibrary.tpre2snowf_cut.restype = ctypes.c_double
        tprecip2snowf_cut = fortlibrary.tpre2snowf_cut(ctypes.byref(c_tprec),
                                                       ctypes.byref(c_tair),
                                                       ctypes.byref(c_tmelt0))
    return tprecip2snowf_cut


def tpre2snowf_lin(tprec, tair, thigh, tlow,
                    library='../src/CISSEMBEL_CBind4Tests.so'):
    """
    C-binding to Fortran function tpre2snowf_lin in mod_physic in CISSEMBEL
    Snowfall rate of total precipitation based on temperature range.

    Parameters
    ----------
    tprec : float or numpy.ndarray
        Total precipitation.
    tair : float or numpy.ndarray
        Air temperature.
    thigh : float
        High temperature threshold above precipitation falls purely as rain.
    tlow : float
        Low temperature threshold above precipitation falls purely as snow.
    library : str or ctypes.CDLL
        Reference to imported shared library containing the C-bindings.

    Returns
    -------
    tprecip2snowf_lin : float or numpy.ndarray
        Snowfall rate from total precipitation.
    """
    if isinstance(library, ctypes.CDLL):
        fortlibrary = library
    else:
        fortlibrary = ctypes.cdll.LoadLibrary(library)

    if isinstance(tprec, (float, int)) and False:
        inum = 0
        # Single c-value
        c_tprec = ctypes.c_double(tprec)
        c_tair = ctypes.c_double(tair)
    else:
        inum = len(tprec)
        # Number conversion, Fortran-order
        c_tprec = numpy.array(tprec, order="F", dtype=float)
        c_tair = numpy.array(tair, order="F", dtype=float)
        # c-pointer
        c_tprec_ptr = c_tprec.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        c_tair_ptr = c_tair.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    c_thigh = ctypes.c_double(thigh)
    c_tlow = ctypes.c_double(tlow)

    if inum:
        # array
        tprecip2snowf_lin = numpy.zeros((inum), order="F", dtype=float)
        tprecip2snowf_lin_ptr = \
            tprecip2snowf_lin.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        fortlibrary.c_tpre2snowf_lin.argtypes = [ctypes.c_int,
                                                 ctypes.POINTER(ctypes.c_double),
                                                 ctypes.POINTER(ctypes.c_double),
                                                 ctypes.c_double,
                                                 ctypes.c_double,
                                                 ctypes.POINTER(ctypes.c_double)]
        fortlibrary.c_tpre2snowf_lin(ctypes.c_int(inum),
                                     c_tprec_ptr,
                                     c_tair_ptr,
                                     c_thigh,
                                     c_tlow,
                                     tprecip2snowf_lin_ptr)
    else:
        # skalar
        fortlibrary.tpre2snowf_lin.restype = ctypes.c_double
        tprecip2snowf_lin = fortlibrary.tpre2snowf_lin(ctypes.byref(c_tprec),
                                                       ctypes.byref(c_tair),
                                                       ctypes.byref(c_thigh),
                                                       ctypes.byref(c_tlow))
    return tprecip2snowf_lin


def hctemp(temp, dz, lapse_rate, library='../src/CISSEMBEL_CBind4Tests.so'):
    """
    C-binding to Fortran function hctemp in mod_physic in CISSEMBEL
    Height correction of air temperature by lapse rate.

    Parameters
    ----------
    temp : float or numpy.ndarray
        Temperature.
    dz : float or numpy.ndarray
        Height difference.
    lapse_rate : float
        lapse rate in atmosphere.
    library : str or ctypes.CDLL
        Reference to imported shared library containing the C-bindings.

    Returns
    -------
    hc_temp : float or numpy.ndarray
        Height-difference corrected temperature.

    """
    if isinstance(library, ctypes.CDLL):
        fortlibrary = library
    else:
        fortlibrary = ctypes.cdll.LoadLibrary(library)

    if isinstance(temp, (float, int)) and False:
        inum = 0
        # Single c-value
        c_temp = ctypes.c_double(temp)
        c_dz = ctypes.c_double(dz)
    else:
        inum = len(temp)
        # Number conversion, Fortran-order
        c_temp = numpy.array(temp, order="F", dtype=float)
        c_dz = numpy.array(dz, order="F", dtype=float)
        # c-pointer
        c_temp_ptr = c_temp.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        c_dz_ptr = c_dz.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    c_lapse_rate = ctypes.c_double(lapse_rate)

    if inum:
        # array
        hc_temp = numpy.zeros((inum), order="F", dtype=float)
        hc_temp_ptr = hc_temp.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        fortlibrary.c_hctemp.argtypes = [ctypes.c_int,
                                         ctypes.POINTER(ctypes.c_double),
                                         ctypes.POINTER(ctypes.c_double),
                                         ctypes.c_double,
                                         ctypes.POINTER(ctypes.c_double)]
        fortlibrary.c_hctemp(ctypes.c_int(inum),
                             c_temp_ptr,
                             c_dz_ptr,
                             c_lapse_rate,
                             hc_temp_ptr)
    else:
        # skalar
        fortlibrary.hctemp.restype = ctypes.c_double
        hc_temp = fortlibrary.hctemp(ctypes.byref(c_temp),
                                     ctypes.byref(c_dz),
                                     ctypes.byref(c_lapse_rate))
    return hc_temp


def hclongwave(lwrad, dz, library='../src/CISSEMBEL_CBind4Tests.so'):
    """
    C-binding to Fortran function hclongwave in mod_physic in CISSEMBEL
    Height correction of longwave/thermal radiation.

    Parameters
    ----------
    lwrad : TYPE
        Longwave / thermal radiation.
    dz : float or numpy.ndarray
        Height difference.
    library : str or ctypes.CDLL
        Reference to imported shared library containing the C-bindings.

    Returns
    -------
    hc_longwave : float or numpy.ndarray
        Height-difference corrected longwave / thermal radiation.

    """
    if isinstance(library, ctypes.CDLL):
        fortlibrary = library
    else:
        fortlibrary = ctypes.cdll.LoadLibrary(library)

    if isinstance(lwrad, float) and False: # TODO todo : add skalar-function to mod_physic
        inum = 0
        # Single c-value
        c_lwrad = ctypes.c_double(lwrad)
        c_dz = ctypes.c_double(dz)
    else:
        inum = len(lwrad)
        # Number conversion, Fortran-order
        c_lwrad = numpy.array(lwrad, order="F", dtype=float)
        c_dz = numpy.array(dz, order="F", dtype=float)
        # c-pointer
        c_lwrad_ptr = c_lwrad.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        c_dz_ptr = c_dz.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

    if inum:
        # array
        hc_longwave = numpy.zeros((inum), order="F", dtype=float)
        hc_longwave_ptr = hc_longwave.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        fortlibrary.c_hclongwave.argtypes = [ctypes.c_int,
                                             ctypes.POINTER(ctypes.c_double),
                                             ctypes.POINTER(ctypes.c_double),
                                             ctypes.POINTER(ctypes.c_double)]
        fortlibrary.c_hclongwave(ctypes.c_int(inum),
                                 c_lwrad_ptr,
                                 c_dz_ptr,
                                 hc_longwave_ptr)
    else:
        # skalar
        fortlibrary.hclongwave.restype = ctypes.c_double
        hc_longwave = fortlibrary.hclongwave(ctypes.byref(c_lwrad),
                                             ctypes.byref(c_dz))
    return hc_longwave


def hcpressure(press, zs_ref, zs_target, library='../src/CISSEMBEL_CBind4Tests.so'):
    """
    C-binding to Fortran function hcpressure in mod_physic in CISSEMBEL
    Height correction of air pressure.

    Parameters
    ----------
    press : float or numpy.ndarray
        air pressure.
    zs_ref : float or numpy.ndarray
        reference topography/orography.
    zs_target : float or numpy.ndarray
        target topography/orography.
    library : str or ctypes.CDLL
        Reference to imported shared library containing the C-bindings.

    Returns
    -------
    hc_pressure : float or numpy.ndarray
        Height-difference corrected air pressure.

    """
    if isinstance(library, ctypes.CDLL):
        fortlibrary = library
    else:
        fortlibrary = ctypes.cdll.LoadLibrary(library)

    if isinstance(press, float) and False: # TODO todo : add skalar-function to mod_physic
        inum = 0
        # Single c-value
        c_press = ctypes.c_double(press)
        c_zs_ref = ctypes.c_double(zs_ref)
        c_zs_target = ctypes.c_double(zs_target)
    else:
        inum = len(press)
        # Number conversion, Fortran-order
        c_press = numpy.array(press, order="F", dtype=float)
        c_zs_ref = numpy.array(zs_ref, order="F", dtype=float)
        c_zs_target = numpy.array(zs_target, order="F", dtype=float)
        # c-pointer
        c_press_ptr = c_press.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        c_zs_ref_ptr = c_zs_ref.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        c_zs_target_ptr = c_zs_target.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

    if inum:
        # array
        if len(press.squeeze().shape) == 2:
            inum = press.shape[0]
            jnum = press.shape[1]
            hc_pressure = numpy.zeros((inum, jnum), order="F", dtype=float)
            hc_pressure_ptr = hc_pressure.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
            fortlibrary.c_hcpressure2D.argtypes = [ctypes.c_int,
                                                   ctypes.c_int,
                                                   ctypes.POINTER(ctypes.c_double),
                                                   ctypes.POINTER(ctypes.c_double),
                                                   ctypes.POINTER(ctypes.c_double),
                                                   ctypes.POINTER(ctypes.c_double)]
            fortlibrary.c_hcpressure2D(ctypes.c_int(inum),
                                       ctypes.c_int(jnum),
                                       c_press_ptr,
                                       c_zs_ref_ptr,
                                       c_zs_target_ptr,
                                       hc_pressure_ptr)
        else:
            hc_pressure = numpy.zeros((inum), order="F", dtype=float)
            hc_pressure_ptr = hc_pressure.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
            fortlibrary.c_hcpressure.argtypes = [ctypes.c_int,
                                                 ctypes.POINTER(ctypes.c_double),
                                                 ctypes.POINTER(ctypes.c_double),
                                                 ctypes.POINTER(ctypes.c_double),
                                                 ctypes.POINTER(ctypes.c_double)]
            fortlibrary.c_hcpressure(ctypes.c_int(inum),
                                     c_press_ptr,
                                     c_zs_ref_ptr,
                                     c_zs_target_ptr,
                                     hc_pressure_ptr)
    else:
        # skalar
        fortlibrary.hcpressure.restype = ctypes.c_double
        hc_pressure = fortlibrary.hcpressure(ctypes.byref(c_press),
                                             ctypes.byref(c_zs_ref),
                                             ctypes.byref(c_zs_target))
    return hc_pressure


def hcprecip_high_desert(precip, zs_ref, zs_target,
                         library='../src/CISSEMBEL_CBind4Tests.so'):
    """
    C-binding to Fortran function hcprecip_high_desert in mod_physic in CISSEMBEL
    Height correction of total precipitation via deserfication effectypes.

    Parameters
    ----------
    precip : TYPE
        DESCRIPTION.
    zs_ref : float or numpy.ndarray
        reference topography/orography.
    zs_target : float or numpy.ndarray
        target topography/orography.
    library : str or ctypes.CDLL
        Reference to imported shared library containing the C-bindings.

    Returns
    -------
    hc_precip_high_desert : TYPE
        Height-difference corrected precipitation due to deserfication process.

    """
    if isinstance(library, ctypes.CDLL):
        fortlibrary = library
    else:
        fortlibrary = ctypes.cdll.LoadLibrary(library)

    if isinstance(precip, float):
        inum = 0
        # Single c-value
        c_precip = ctypes.c_double(precip)
        c_zs_ref = ctypes.c_double(zs_ref)
        c_zs_target = ctypes.c_double(zs_target)
    else:
        inum = len(precip)
        # Number conversion, Fortran-order
        c_precip = numpy.array(precip, order="F", dtype=float)
        c_zs_ref = numpy.array(zs_ref, order="F", dtype=float)
        c_zs_target = numpy.array(zs_target, order="F", dtype=float)
        # c-pointer
        c_precip_ptr = c_precip.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        c_zs_ref_ptr = c_zs_ref.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        c_zs_target_ptr = c_zs_target.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

    if inum:
        # array
        if len(precip.squeeze().shape) == 2:
            inum = precip.shape[0]
            jnum = precip.shape[1]
            hc_precip_high_desert = numpy.zeros((inum, jnum), order="F", dtype=float)
            hc_precip_high_desert_ptr = \
                hc_precip_high_desert.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
            fortlibrary.c_hcprecip_high_desert2D.argtypes = [ctypes.c_int,
                                                             ctypes.c_int,
                                                             ctypes.POINTER(ctypes.c_double),
                                                             ctypes.POINTER(ctypes.c_double),
                                                             ctypes.POINTER(ctypes.c_double),
                                                             ctypes.POINTER(ctypes.c_double)]
            fortlibrary.c_hcprecip_high_desert2D(ctypes.c_int(inum),
                                                 ctypes.c_int(jnum),
                                                 c_precip_ptr,
                                                 c_zs_ref_ptr,
                                                 c_zs_target_ptr,
                                                 hc_precip_high_desert_ptr)
        else:
            hc_precip_high_desert = numpy.zeros((inum), order="F", dtype=float)
            hc_precip_high_desert_ptr = \
                hc_precip_high_desert.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
            fortlibrary.c_hcprecip_high_desert.argtypes = [ctypes.c_int,
                                                           ctypes.POINTER(ctypes.c_double),
                                                           ctypes.POINTER(ctypes.c_double),
                                                           ctypes.POINTER(ctypes.c_double),
                                                           ctypes.POINTER(ctypes.c_double)]
            fortlibrary.c_hcprecip_high_desert(ctypes.c_int(inum),
                                               c_precip_ptr,
                                               c_zs_ref_ptr,
                                               c_zs_target_ptr,
                                               hc_precip_high_desert_ptr)
    else:
        # skalar
        fortlibrary.hcprecip_high_desert.restype = ctypes.c_double
        hc_precip_high_desert = \
            fortlibrary.hcprecip_high_desert(ctypes.byref(c_precip),
                                             ctypes.byref(c_zs_ref),
                                             ctypes.byref(c_zs_target))
    return hc_precip_high_desert


def snowlayers2pressure(thickness, density, patm,
                        library='../src/CISSEMBEL_CBind4Tests.so'):
    """
    C-binding to Fortran function/subroutine snowlayers2pressure in mod_physic in CISSEMBEL
    Pressure at snow layer(s).

    Parameters
    ----------
    thickness : float or numpy.ndarray
        layer thickness(es).
    density : float or numpy.ndarray
        Density of layer(s).
    patm : float
        Atmospheric standard pressure.
    library : str or ctypes.CDLL
        Reference to imported shared library containing the C-bindings.

    Returns
    -------
    snowlayers2pressure : float or numpy.ndarray
        Pressure at the layers center.

    """
    if isinstance(library, ctypes.CDLL):
        fortlibrary = library
    else:
        fortlibrary = ctypes.cdll.LoadLibrary(library)

    if isinstance(thickness, (float, int)):
        inum = 0
        # Single c-value
        c_thickness = ctypes.c_double(thickness)
        c_density = ctypes.c_double(density)
    else:
        inum = len(thickness)
        # Number conversion, Fortran-order
        c_thickness = numpy.array(thickness, order="F", dtype=float)
        c_density = numpy.array(density, order="F", dtype=float)
        # c-pointer
        c_thickness_ptr = c_thickness.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        c_density_ptr = c_density.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    c_patm = ctypes.c_double(patm)

    if inum:
        # array
        snowlay2press = numpy.zeros((inum), order="F", dtype=float)
        snowlay2press_ptr = \
            snowlay2press.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        fortlibrary.c_snowlayers2pressure.argtypes = [ctypes.c_int,
                                                      ctypes.POINTER(ctypes.c_double),
                                                      ctypes.POINTER(ctypes.c_double),
                                                      ctypes.POINTER(ctypes.c_double),
                                                      ctypes.c_double]
        fortlibrary.c_snowlayers2pressure(ctypes.c_int(inum),
                                          c_thickness_ptr,
                                          c_density_ptr,
                                          snowlay2press_ptr,
                                          c_patm)
    else:
        # skalar
        fortlibrary.snowlayers2pressure.restype = ctypes.c_double
        snowlay2press = fortlibrary.snowlayers2pressure(ctypes.byref(c_thickness),
                                                              ctypes.byref(c_density),
                                                              ctypes.byref(c_patm))
    return snowlay2press


def densification_pressure(dtime, pressure, temperature, density, scont, icont,
                           wcont, Tfreezing, library='../src/CISSEMBEL_CBind4Tests.so'):
    """
    C-binding to Fortran function densification_pressure in mod_physic in CISSEMBEL
    Density evolution of snow layer(s).

    Parameters
    ----------
    dtime : float
        Time step.
    pressure : float or numpy.ndarray
        Pressure at snow layer(s).
    temperature : float or numpy.ndarray
        Temperature of snow layer(s).
    density : float or numpy.ndarray
        Density of snow layer(s).
    scont : float or numpy.ndarray
        Snow content of snow layer(s).
    icont : float or numpy.ndarray
        Ice content of snow layer(s).
    wcont : float or numpy.ndarray
        Water content of snow layer(s).
    Tfreezing : float or numpy.ndarray
        Freezing point temperature content of snow layer(s).
    library : str or ctypes.CDLL
        Reference to imported shared library containing the C-bindings.

    Returns
    -------
    densification_pres : float or numpy.ndarray
        Density of snow layer(s) .

    """
    if isinstance(library, ctypes.CDLL):
        fortlibrary = library
    else:
        fortlibrary = ctypes.cdll.LoadLibrary(library)

    if isinstance(pressure, (float, int)):
        inum = 0
        # Single c-value
        c_pressure = ctypes.c_double(pressure)
        c_temperature = ctypes.c_double(temperature)
        c_density = ctypes.c_double(density)
        c_scont = ctypes.c_double(scont)
        c_icont = ctypes.c_double(icont)
        c_wcont = ctypes.c_double(wcont)
        c_Tfreezing = ctypes.c_double(Tfreezing)
    else:
        inum = len(pressure)
        # Number conversion, Fortran-order
        c_pressure = numpy.array(pressure, order="F", dtype=float)
        c_temperature = numpy.array(temperature, order="F", dtype=float)
        c_density = numpy.array(density, order="F", dtype=float)
        c_scont = numpy.array(scont, order="F", dtype=float)
        c_icont = numpy.array(icont, order="F", dtype=float)
        c_wcont = numpy.array(wcont, order="F", dtype=float)
        c_Tfreezing = numpy.array(Tfreezing, order="F", dtype=float)
        # c-pointer
        c_pressure_ptr = c_pressure.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        c_temperature_ptr = c_temperature.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        c_density_ptr = c_density.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        c_scont_ptr = c_scont.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        c_icont_ptr = c_icont.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        c_wcont_ptr = c_wcont.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        c_Tfreezing_ptr = c_Tfreezing.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    c_dtime = ctypes.c_double(dtime)

    if inum:
        # array
        densification_pres = numpy.zeros((inum), order="F", dtype=float)
        densification_pres_ptr = \
            densification_pres.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        fortlibrary.c_densification_pressure.argtypes = [ctypes.c_int,
                                                         ctypes.c_double,
                                                         ctypes.POINTER(ctypes.c_double),
                                                         ctypes.POINTER(ctypes.c_double),
                                                         ctypes.POINTER(ctypes.c_double),
                                                         ctypes.POINTER(ctypes.c_double),
                                                         ctypes.POINTER(ctypes.c_double),
                                                         ctypes.POINTER(ctypes.c_double),
                                                         ctypes.POINTER(ctypes.c_double),
                                                         ctypes.POINTER(ctypes.c_double)]
        fortlibrary.c_densification_pressure(ctypes.c_int(inum),
                                             c_dtime,
                                             c_pressure_ptr,
                                             c_temperature_ptr,
                                             c_density_ptr,
                                             c_scont_ptr,
                                             c_icont_ptr,
                                             c_wcont_ptr,
                                             c_Tfreezing_ptr,
                                             densification_pres_ptr)
    else:
        # skalar
        fortlibrary.densification_pressure.restype = ctypes.c_double
        densification_pres = fortlibrary.densification_pressure(ctypes.byref(c_dtime),
                                                                ctypes.byref(c_pressure),
                                                                ctypes.byref(c_temperature),
                                                                ctypes.byref(c_density),
                                                                ctypes.byref(scont),
                                                                ctypes.byref(icont),
                                                                ctypes.byref(wcont),
                                                                ctypes.byref(Tfreezing))
    return densification_pres


def waterdepth2pressure(waterdepth, rho_water, patm,
                        library='../src/CISSEMBEL_CBind4Tests.so'):
    """
    C-binding to Fortran function waterdepth2pressure in mod_physic in CISSEMBEL
    Pressure at given water depth.

    Parameters
    ----------
    waterdepth : float or numpy.ndarray
        Water depth.
    rho_water : float or numpy.ndarray
        Density of water at given depth(s).
    patm : float
        Atmospheric standard pressure.
    library : str or ctypes.CDLL
        Reference to imported shared library containing the C-bindings.

    Returns
    -------
    water_column2pressure : float or numpy.ndarray
        Pressure at given water depth(s).

    """
    if isinstance(library, ctypes.CDLL):
        fortlibrary = library
    else:
        fortlibrary = ctypes.cdll.LoadLibrary(library)

    if isinstance(waterdepth, (float, int)):
        inum = 0
        # Single c-value
        c_waterdepth = ctypes.c_double(waterdepth)
        c_rho_water = ctypes.c_double(rho_water)
    else:
        inum = len(waterdepth)
        # Number conversion, Fortran-order
        c_waterdepth = numpy.array(waterdepth, order="F", dtype=float)
        c_rho_water = numpy.array(rho_water, order="F", dtype=float)
        # c-pointer
        c_waterdepth_ptr = c_waterdepth.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        c_rho_water_ptr = c_rho_water.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    c_patm = ctypes.c_double(patm)

    if inum:
        # array
        water_column2pressure = numpy.zeros((inum), order="F", dtype=float)
        water_column2pressure_ptr = \
            water_column2pressure.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        fortlibrary.c_waterdepth2pressure.argtypes = [ctypes.c_int,
                                                  ctypes.POINTER(ctypes.c_double),
                                                  ctypes.POINTER(ctypes.c_double),
                                                  ctypes.c_double,
                                                  ctypes.POINTER(ctypes.c_double)]
        fortlibrary.c_waterdepth2pressure(ctypes.c_int(inum),
                                          c_waterdepth_ptr,
                                          c_rho_water_ptr,
                                          c_patm,
                                          water_column2pressure_ptr)
    else:
        # skalar
        fortlibrary.waterdepth2pressure.restype = ctypes.c_double
        water_column2pressure = fortlibrary.waterdepth2pressure(ctypes.byref(c_waterdepth),
                                                              ctypes.byref(c_rho_water),
                                                              ctypes.byref(c_patm))
    return water_column2pressure


def height2pressure(height, rho, patm, library='../src/CISSEMBEL_CBind4Tests.so'):
    """
    C-binding to Fortran function height2pressure in mod_physic in CISSEMBEL
    Pressure below a column of matter with density rho.

    Parameters
    ----------
    height : float or numpy.ndarray
        Height of mass column of density rho.
    rho : float or numpy.ndarray
        Density of air.
    patm : float
        Atmospheric standard pressure.
    library : str or ctypes.CDLL
        Reference to imported shared library containing the C-bindings.

    Returns
    -------
    mass_column2pressure : float or numpy.ndarray
        Pressure below column with height(s).

    """
    if isinstance(library, ctypes.CDLL):
        fortlibrary = library
    else:
        fortlibrary = ctypes.cdll.LoadLibrary(library)

    if isinstance(height, (float, int)):
        inum = 0
        # Single c-value
        c_height = ctypes.c_double(height)
        c_rho = ctypes.c_double(rho)
    else:
        inum = len(height)
        # Number conversion, Fortran-order
        c_height = numpy.array(height, order="F", dtype=float)
        c_rho = numpy.array(rho, order="F", dtype=float)
        # c-pointer
        c_height_ptr = c_height.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        c_rho_ptr = c_rho.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    c_patm = ctypes.c_double(patm)

    if inum:
        # array
        mass_column2pressure = numpy.zeros((inum), order="F", dtype=float)
        mass_column2pressure_ptr = \
            mass_column2pressure.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        fortlibrary.c_height2pressure.argtypes = [ctypes.c_int,
                                                  ctypes.POINTER(ctypes.c_double),
                                                  ctypes.POINTER(ctypes.c_double),
                                                  ctypes.c_double,
                                                  ctypes.POINTER(ctypes.c_double)]
        fortlibrary.c_height2pressure(ctypes.c_int(inum),
                                      c_height_ptr,
                                      c_rho_ptr,
                                      c_patm,
                                      mass_column2pressure_ptr)
    else:
        # skalar
        fortlibrary.height2pressure.restype = ctypes.c_double
        mass_column2pressure = fortlibrary.height2pressure(ctypes.byref(c_height),
                                                      ctypes.byref(c_rho),
                                                      ctypes.byref(c_patm))
    return mass_column2pressure


def TfrezSeaWater(pres, salt, library='../src/CISSEMBEL_CBind4Tests.so'):
    """
    C-binding to Fortran function TfrezSeaWater in mod_physic in CISSEMBEL
    Freezing point temperature of sea water.

    Parameters
    ----------
    pres : float or numpy.ndarray
        Pressure in water.
    salt : float or numpy.ndarray
        Salinity.
    library : str or ctypes.CDLL
        Reference to imported shared library containing the C-bindings.

    Returns
    -------
    tfrez_seawater : float or numpy.ndarray
        Freezing point temperature.

    """
    if isinstance(library, ctypes.CDLL):
        fortlibrary = library
    else:
        fortlibrary = ctypes.cdll.LoadLibrary(library)

    if isinstance(pres, (float, int)):
        inum = 0
        # Single c-value
        c_pres = ctypes.c_double(pres)
        c_salt = ctypes.c_double(salt)
    else:
        inum = len(pres)
        # Number conversion, Fortran-order
        c_pres = numpy.array(pres, order="F", dtype=float)
        c_salt = numpy.array(salt, order="F", dtype=float)
        # c-pointer
        c_pres_ptr = c_pres.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        c_salt_ptr = c_salt.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

    if inum:
        # array
        tfrez_seawater = numpy.zeros((inum), order="F", dtype=float)
        tfrez_seawater_ptr = tfrez_seawater.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        fortlibrary.c_TfrezSeaWater.argtypes = [ctypes.c_int,
                                                ctypes.POINTER(ctypes.c_double),
                                                ctypes.POINTER(ctypes.c_double),
                                                ctypes.POINTER(ctypes.c_double)]

        fortlibrary.c_TfrezSeaWater(ctypes.c_int(inum),
                                    c_pres_ptr,
                                    c_salt_ptr,
                                    tfrez_seawater_ptr)
    else:
        # skalar
        fortlibrary.TfrezSeaWater.restype = ctypes.c_double
        tfrez_seawater = fortlibrary.TfrezSeaWater(ctypes.byref(c_pres),
                                                   ctypes.byref(c_salt))
    return tfrez_seawater


def TfrezSoilWater(salt, library='../src/CISSEMBEL_CBind4Tests.so'):
    """
    C-binding to Fortran function TfrezSoilWater in mod_physic in CISSEMBEL
    Freezing point temperature of soil water.

    Parameters
    ----------
    salt : float or numpy.ndarray
        Salinity.
    library : str or ctypes.CDLL
        Reference to imported shared library containing the C-bindings.

    Returns
    -------
    tfrez_soilwater : float or numpy.ndarray
        Freezing point temperature.

    """
    if isinstance(library, ctypes.CDLL):
        fortlibrary = library
    else:
        fortlibrary = ctypes.cdll.LoadLibrary(library)

    if isinstance(salt, (float, int)):
        inum = 0
        # Single c-value
        c_salt = ctypes.c_double(salt)
    else:
        inum = len(salt)
        # Number conversion, Fortran-order
        c_salt = numpy.array(salt, order="F", dtype=float)
        # c-pointer
        c_salt_ptr = c_salt.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

    if inum:
        # array
        tfrez_soilwater = numpy.zeros((inum), order="F", dtype=float)
        tfrez_soilwater_ptr = \
            tfrez_soilwater.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        fortlibrary.c_TfrezSoilWater.argtypes = [ctypes.c_int,
                                                 ctypes.POINTER(ctypes.c_double),
                                                 ctypes.POINTER(ctypes.c_double)]
        fortlibrary.c_TfrezSoilWater(ctypes.c_int(inum),
                                     c_salt_ptr,
                                     tfrez_soilwater_ptr)
    else:
        # skalar
        fortlibrary.TfrezSoilWater.restype = ctypes.c_double
        tfrez_soilwater = fortlibrary.TfrezSoilWater(ctypes.byref(c_salt))
    return tfrez_soilwater


def compute_runoff_time(ElevGrad, library='../src/CISSEMBEL_CBind4Tests.so'):
    """
    C-binding to Fortran function compute_runoff_time in mod_physic in CISSEMBEL
    Local runoff time-scale.

    Parameters
    ----------
    ElevGrad : float or numpy.ndarray
        Surface elevation gradient.
    library : str or ctypes.CDLL
        Reference to imported shared library containing the C-bindings.

    Returns
    -------
    runoff_time : float or numpy.ndarray
        Freezing point temperature.

    """
    if isinstance(library, ctypes.CDLL):
        fortlibrary = library
    else:
        fortlibrary = ctypes.cdll.LoadLibrary(library)

    if isinstance(ElevGrad, (float, int)):
        inum = 0
        # Single c-value
        c_ElevGrad = ctypes.c_double(ElevGrad)
    else:
        inum = len(ElevGrad)
        # Number conversion, Fortran-order
        c_ElevGrad = numpy.array(ElevGrad, order="F", dtype=float)
        # c-pointer
        c_ElevGrad_ptr = c_ElevGrad.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

    if inum:
        # array
        runoff_time = numpy.zeros((inum), order="F", dtype=float)
        runoff_time_ptr = \
            runoff_time.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        fortlibrary.c_compute_runoff_time.argtypes = [ctypes.c_int,
                                                 ctypes.POINTER(ctypes.c_double),
                                                 ctypes.POINTER(ctypes.c_double)]
        fortlibrary.c_compute_runoff_time(ctypes.c_int(inum),
                                     c_ElevGrad_ptr,
                                     runoff_time_ptr)
    else:
        # skalar
        fortlibrary.compute_runoff_time.restype = ctypes.c_double
        runoff_time = fortlibrary.compute_runoff_time(ctypes.byref(c_ElevGrad))
    return runoff_time


# -----------------------------------------------------------
#
# Main program
#
if __name__ == '__main__':
    import ctypes
    import numpy
    import values_mod_param as physical

    #
    # Example code
    #
    print('Load physical constant')
    physical.constant()
    RHO_ICE = physical.constant.rho_ice
    RHO_SNOW = physical.constant.rho_snow
    CDENSITY = physical.constant.Cdens

    #
    # Load the shared library
    #
    LIBRARY_NAME = '../src/CISSEMBEL_CBind4Tests.so' # INCLUDING path, e.g., ./

    print('= Load shared library containing C-bindings of Fortran code "'
          +LIBRARY_NAME+'"')
    FORTLIBRARY = ctypes.cdll.LoadLibrary(LIBRARY_NAME)
else:
    import ctypes
    import numpy
