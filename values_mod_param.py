#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  3 14:58:40 2021

-- Some variables and parameters of mod_param.F08

@author: Christian Rodehacke, DMI Copenhagen, Denmark
"""

def add_this_arg(func):
    """
    Functions are objects in Python and can have arbitrary attributes assigned to them.

    A generic decorator that adds a this argument to each call to the decorated
    function. This additional argument will give functions a way to reference
    themselves without needing to explicitly embed (hardcode) their name into
    the rest of the definition and is similar to the instance argument that
    class methods automatically receive as their first argument which is
    usually named `self`. You may a different one to avoid confusion, such as
    `this`.

    Parameters
    ----------
    func : function
        function call.

    Returns
    -------
    function
        decorated function call.

    """
    def wrapped(*args, **kwargs):
        return func(wrapped, *args, **kwargs)
    return wrapped

@add_this_arg
def constant(this=None):
    '''
    Define various physical constants

    Parameters
    ----------
    this : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    # Physical constant
    this.gravitational_acceleration = 9.80665 # Gravity constant (m s-2)
    this.stefan_boltzmann_constant = 5.670374419e-8 # Stefan Boltzman constant (W m-2 K-4)
    this.gas_constant = 8.31446261815324 # gas constant (J K-1 mol-1)

    # Heat capacity
    this.cp_water = 4186.0    # Heat capacity of fresh water (J kg-1 K-1)
    this.cp_snow  = 2000.0    # Heat capacity of snow (J kg-1 K-1)
    this.cp_ice   = 2090.0    # Heat capacity of ice (J kg-1 K-1)

    # Thermal conductivity
    this.cond_fwater= 0.56789 # Thermal conductivity of water, (W m-1 K-1)
    this.cond_ice = 2.4567    # Thermal conductivity of ice, (W m-1 K-1)

    # Density
    this.rho_fwater  = 1000.0 # Density of fresh water (kg m-3)
    this.rho_snow    =  300.0 # Density of fresh snow (kg m-3)
    this.rho_ice     =  910.0 # Density of ice/glaciers (kg m-3)
    this.rho_ground  = 1748.0 # Density of the ground (kg m-3)
    this.rho_seawater= 1028.0 # Density of sea water (kg m-3)

    # Density profile
    this.Cdens = 0.024

    # Atmospheric
    this.patm_surf   = 1000.0e2 # Atmospheric surface pressure (Pa)
    this.lapse_rate  = -6.5e-3  # Atmospheric lapse rate (K m-1)
    this.gradlw  = -75.0/2600.0 # Vertical longwave gradient (W m-2), Marty (2002)

    # Fresh water
    this.Tmelt_fw = 273.15      # melting/freezing temperature of freshwater (K)

# -----------------------------------------------------------
#
# Main program
#
if __name__ == '__main__':
    constant()

    # # Physical constant
    # gravitational_acceleration = 9.80665 # Gravity constant (m s-2)
    # stefan_boltzmann_constant = 5.670374419e-8 # Stefan Boltzman constant (W m-2 K-4)
    # gas_constant = 8.31446261815324 # gas constant (J K-1 mol-1)

    # # Heat capacity
    # cp_water = 4186.0    # Heat capacity of fresh water (J kg-1 K-1)
    # cp_snow  = 2000.0    # Heat capacity of snow (J kg-1 K-1)
    # cp_ice   = 2090.0    # Heat capacity of ice (J kg-1 K-1)

    # # Thermal conductivity
    # cond_fwater= 0.56789 # Thermal conductivity of water, (W m-1 K-1)
    # cond_ice = 2.4567    # Thermal conductivity of ice, (W m-1 K-1)

    # # Density
    # rho_fwater  = 1000.0 # Density of fresh water (kg m-3)
    # rho_snow    =  300.0 # Density of fresh snow (kg m-3)
    # rho_ice     =  910.0 # Density of ice/glaciers (kg m-3)
    # rho_ground  = 1748.0 # Density of the ground (kg m-3)
    # rho_seawater= 1028.0 # Density of sea water (kg m-3)

    # # Density profile
    # Cdens = 0.024
