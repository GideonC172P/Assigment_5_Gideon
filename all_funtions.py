# -*- coding: utf-8 -*-
"""
Created on Fri Jul 17 18:28:30 2020

@author: Gideon
"""
def albedo(alphaw, alphab, alphag, aw, ab, ag):
    """Caluclate the average planetary albedo. 
    alphaw*aw + alphab*ab + alphag*ag
    
    Arguments
    ---------
    alphaw : float
        Surface fraction of white daisies.
    alphab : float
        Surface fraction of black daisies.
    alphag : float
        Surface fraction of bare ground.
    aw = 0.7 : float
        Albedo of white daisies.
    ab = 0.8 : float
        Albedo of black dasies.
    ag = 0.5 : float
        Albedo of bare ground
    """
    return alphaw*aw + alphab*ab + alphag*ag


def daisy_replicator(alpha, alphag, beta, gamma):
    """Calculate the rate of change of a replicator (daisies). The formula is:
    a*[ag*b-g]
    which describes a logistic growth and constant death within a limited 
    resource system over time t.
    
    Arguments
    ---------
    alpha :  float
        alpha -- the fraction of surface covered by daisies.
    alphag : float
        alpha_g -- The fraction of available bare ground.
    beta :  float
        beta -- the birth rate.
    gamma : float
        gamma -- the death rate.
    """
    return alpha*(alphag*beta-gamma)


def beta(temp, optimum=295.65, tmin=278.15, tmax=313.15, k=0.003265):
    """Calculate the environment parameter beta. This is like a birth rate
    that reaches a maximum at an optimal temperature, and a range of birth
    rates is specified by a parabolic width parameter of k=0.003265.
    1 - k*(temp-optimum)**2
    
    Arguments
    ---------
    temp : float
        The environment temperature experienced by the daisies
    
    Keyword arguments
    -----------------
    optimum = 295.65 : float
        The optimum temperature for daisy birthrate.
    tmin = 278.15 : float
        Maximum temperature for birthrate.
    tmax = 313.15 : float
        Minimum temperature for birthrate.
    k = 0.003265 : float
        Quadratic width parameter.
    """
    if (temp>tmin)&(temp<tmax):
        return 1 - k*(temp-optimum)**2
    else:
        return 0


def euler(initial, tendency, h=1):
    """Integrate forward in time using Euler's method of numerical integration.
    initial + h*tendency
    
    Arguments
    ---------
    initial :  float
        The initial state.
    tendency : float
        The rate of change in the initial state.
    
    Keyword arguments
    -----------------
    h = 1 : float
        The timestep duration.
    """
    return initial + h*tendency


def local_temp(A, albedo, T, q=30):
    """Calculate local temperature experienced by a particular daisy type or 
    ground cover. This is a simplified version of the original.
    q*(A-albedo)+T
    
    Arguments
    ---------
    A : float
        Planetary albedo
    alpha : float
        Albedo of daisy type
    T : float
        Planetary temperature.
    
    Keyword Arguments
    -----------------
    q = 30 : float
        Heat transfer coefficient
    """
    return q*(A-albedo)+T


def planetary_temp(S, A, L=1.0):
    """Calculate the planetary temperature.
    SL(1-A) = sT**4
    
    Arguments
    ---------
    S : float
        Incident solar energy.
    A : float
        Planetary albedo.

    Keyword Arguments
    -----------------
    L = 1.0 : float
        Normalised stellar luminosity.
    """
    sigma = 5.67032e-8 # Stephan-Bolzmann constant.
    return ((S*L*(1-A))/sigma)**(1/4.)


