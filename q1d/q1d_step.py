#!/usr/bin/env python

from __future__ import division
from math import *
import numpy as np
from numpy.fft import fft, fftshift, fftfreq

def fermi(x, smearing):
    """Return Fermi function"""
    return 1./(1. + np.exp(x/smearing))

def step(x, delta, smearing):
    """Return smoothed step-function Fermi(x-Delta)-Fermi(x)"""
    return fermi(x-delta,smearing) - fermi(x,smearing)

def build_wire(x, heights, delta, smearing):
    """Return sum_(n=-N,N) alpha_n step(x-n*Delta)"""
    N_module = heights.size
    wire = np.zeros(x.size)
    for n in np.arange(N_module):
        wire = ( wire + heights[n]*
                 step(x - (n-N_module/2.)*delta, delta, smearing) )
    return wire

def powerspectrum(data, dx):
    """Return power-spectrum of input signal"""
    powerspec = np.abs(fftshift(fft(data))*dx)**2
    freq = 2.*pi*fftfreq(data.size, dx)
    return freq, powerspec

def AGS(k, heights, delta, smearing):
    """Return roughness-height power spectrum W(k)"""
    N_module = np.arange(heights.size)
    # Take correlation of alphas into account:
    omega = ([ np.exp(-1j*n*k*delta)*heights[n] for n in N_module ])
    omega = np.sum(omega, axis=0)
    omega = np.abs(omega)**2 / heights.size
    return (1./delta * (2.*pi*smearing*
            np.sinh(k*pi*smearing)**(-1)*np.sin(k*delta/2.))**2) * omega

def SGS(k, heights, delta, smearing):
    """Return roughness-height power spectrum S(k)"""
    N_module = np.arange(heights.size)
    # a[n-1] and a[n+1] = 0 for n=N and n=0
    a = np.concatenate([ [0],heights,[0] ])  
    # Take correlation of alphas into account:
    omega = ([ np.exp(-1j*n*k*delta)*(a[n]*(a[n]-a[n+1])*
               np.exp(-1j*k*delta) + a[n]*(a[n]-a[n-1])) for n in N_module ])
    omega = np.sum(omega, axis=0)  
    omega = np.abs(omega)**2 / heights.size
    return ( 1./delta / 72. * (k*pi*(1.+k**2*smearing**2)*
             np.sinh(k*pi*smearing)**(-1))**2 ) * omega

def transmission(n, d, L, sigma, k, heights, delta, smearing):
    """
    Return transmission T based on analytical expressions for W and S
    """
    # only symmetric wire geometry considered yet
    invLbAGS=(4.*sigma**2 / d**6 *(n*pi)**4 /
              k**2 * AGS(2*k,heights,delta,smearing))
    invLbSGS=(0.5*(sigma/d*pi*n)**4 / k**2 *
              (1./3. + 1./(pi*n)**2)**2 * SGS(2*k,heights,delta,smearing))
    return np.exp(-L*(invLbAGS+invLbSGS)), invLbAGS, invLbSGS
