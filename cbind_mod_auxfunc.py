#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 21 20:00:17 CET 2022

C-binding to Fortran functions and subroutines in mod_auxfunc as part of the
Copenhagen Ice Snow Surface Energy and Mass Balance modEL (CISSEMBEL)

@author: Christian Rodehacke, DMI Copenhagen, Denmark
"""
# TODO : Improve python convention, refactor, warning

def rescale_field(reference, field2rescale, mask,
                  library='../src/CISSEMBEL_CBind4Tests.so'):
    """
    C-binding to Fortran function rescale_field in mod_physic in CISSEMBEL
    Factor to rescale a field

    Parameters
    ----------
    reference : numpy.ndarray
        Reference field, source field.
    field2rescale : numpy.ndarray
        DField that shall be rescaled.
    mask : numpy.ndarray of bool/logical
        Mask which points shall be included in the rescaling.
    library : str or ctypes.CDLL
        Reference to imported shared library containing the C-bindings.

    Returns
    -------
    rescale_factor rescaled_field : float
        Rescaling factor.

    """
    if isinstance(library, ctypes.CDLL):
        fortlibrary = library
    else:
        fortlibrary = ctypes.cdll.LoadLibrary(library)

    inum = len(reference)
    if inum <=1:
        print('** cbind_mod_physic:rescale_field'+
              ' Length of profile too short: '+str(inum)+'; need at least 2')
        return -9999999999999.9

    # Number conversion, Fortran-order
    c_reference = numpy.array(reference, order="F", dtype=float)
    c_field2rescale = numpy.array(field2rescale, order="F", dtype=float)
    c_mask = numpy.array(mask, order="F", dtype=bool)
    # c-pointer
    c_reference_ptr = c_reference.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    c_field2rescale_ptr = c_field2rescale.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    c_mask_ptr = c_mask.ctypes.data_as(ctypes.POINTER(ctypes.c_bool))

    # array function to skalar
    if len(reference.squeeze().shape) == 2:
        inum = reference.shape[0]
        jnum = reference.shape[1]
        fortlibrary.rescale_field2D.argtypes = [ctypes.c_int,
                                                ctypes.c_int,
                                                ctypes.POINTER(ctypes.c_double),
                                                ctypes.POINTER(ctypes.c_double),
                                                ctypes.POINTER(ctypes.c_bool)]
        fortlibrary.rescale_field2D.restype = ctypes.c_double
        rescale_factor = fortlibrary.rescale_field2D(ctypes.c_int(inum),
                                                     ctypes.c_int(jnum),
                                                     c_reference_ptr,
                                                     c_field2rescale_ptr,
                                                     c_mask_ptr)
    else:
        fortlibrary.rescale_field.argtypes = [ctypes.c_int,
                                              ctypes.POINTER(ctypes.c_double),
                                              ctypes.POINTER(ctypes.c_double),
                                              ctypes.POINTER(ctypes.c_bool)]
        fortlibrary.rescale_field.restype = ctypes.c_double
        rescale_factor = fortlibrary.rescale_field(ctypes.c_int(inum),
                                                   c_reference_ptr,
                                                   c_field2rescale_ptr,
                                                   c_mask_ptr)
    return rescale_factor


def values_at_depth_level(field, depth, depth_subsurface=10.,
                                 library='../src/CISSEMBEL_CBind4Tests.so'):
    """
    C-binding to Fortran function values_at_depth_level in mod_physic in CISSEMBEL
    Interpolation of property field profile on depth_subsurface.

    Parameters
    ----------
    field : numpy.ndarray
        propertiy profile.
    depth : numpy.ndarray
        depth profile.
    depth_subsurface : float
        Target depth. The default is 10..
    library : str or ctypes.CDLL
        Reference to imported shared library containing the C-bindings.

    Returns
    -------
    value_subsurface : float
        Property value at required depth.

    """
    if isinstance(library, ctypes.CDLL):
        fortlibrary = library
    else:
        fortlibrary = ctypes.cdll.LoadLibrary(library)

    inum = len(field)
    if inum <=1:
        print('** cbind_mod_physic:values_at_depth_level'+
              ' Length of profile too short: '+str(inum)+'; need at least 2')

    # Number conversion, Fortran-order
    c_field = numpy.array(field, order="F", dtype=float)
    c_depth = numpy.array(depth, order="F", dtype=float)
    # c-pointer
    c_field_ptr = c_field.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    c_depth_ptr = c_depth.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    # Single c-value
    c_depth_subsurface = ctypes.c_double(depth_subsurface)
    # array function to skalar
    fortlibrary.values_at_depth_level.argtypes = [ctypes.c_int,
                                                  ctypes.POINTER(ctypes.c_double),
                                                  ctypes.POINTER(ctypes.c_double),
                                                  ctypes.c_double]
    fortlibrary.values_at_depth_level.restype = ctypes.c_double
    value_subsurface = fortlibrary.values_at_depth_level(ctypes.c_int(inum),
                                                         c_field_ptr,
                                                         c_depth_ptr,
                                                         c_depth_subsurface)
    return value_subsurface

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
