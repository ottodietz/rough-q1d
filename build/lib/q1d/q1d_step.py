#!/usr/bin/env python

import sys
import argparse
from math import *
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt


# -----------------------------------------------------------------------------
# Function definitions
# -----------------------------------------------------------------------------

def loadAlphas(alphaFile):
    """ this function is redundant """
    return [ (x - x.mean()) for x in ( (0.8 + 0.2*y)/100 for y in np.loadtxt(alphaFile,unpack=True) ) ]

def smooth_step(x,delta,smearing):
    """
    Return smoothed step-function Fermi(x-Delta)-Fermi(x)
    """
    return 1./(1.+np.exp((x-delta)/smearing))-1./(1.+np.exp(x/smearing))

def build_wire(x,alpha,delta,smearing):
    """
    Return sum_(n=-N,N) alpha_n step(x-n*Delta)
    """
    module_number = alpha.size
    temp_wire = np.zeros(x.size)
    for n in np.arange(module_number):
        temp_wire = (temp_wire + alpha[n]*
                     step(x-(n-module_number/2.)*delta,delta,smearing))
    return temp_wire

def powerspectrum(data):
    """
    Return power-spectrum of input signal
    """
    dx = abs(data[2] - data[1])
    FFT = np.fft.fftshift(abs(np.fft.fft(data)*dx))**2
    freq = np.fft.fftfreq(data.size, dx)
    return freq, FFT

def AGS(k,alpha,delta,smearing):
    """
    Return roughness-height power spectrum W(k)
    """
    module_index = np.arange(alpha.size)
    # Take correlation of alphas into account:
    omega = ([ np.exp(-1j*n*k*delta)*alpha[n] for n in module_index])
    omega = np.sum(omega,axis=0)
    omega = abs(omega)**2 / alpha.size
    return (1./delta * (2.*pi*smearing*
            np.sinh(k*pi*smearing)**(-1)*np.sin(k*delta/2.))**2) * omega

def SGS(k,alpha,delta,smearing):
    """
    Return roughness-height power spectrum S(k)
    """
    module_index = np.arange(alpha.size)
    # a[n-1] and a[n+1] = 0 for n=N and n=0
    a = np.concatenate([ [0],alpha,[0] ])  
    # Take correlation of alphas into account:
    omega = ([ np.exp(-1j*n*k*delta)*(a[n]*(a[n]-a[n+1])*
               np.exp(-1j*k*delta) + a[n]*(a[n]-a[n-1])) for n in module_index])
    omega = np.sum(omega,axis=0)  
    omega = abs(omega)**2 / alpha.size 
    return ( 1./delta / 72. * (k*pi*(1.+k**2*smearing**2)*
             np.sinh(k*pi*smearing)**(-1))**2 ) * omega

def T_analytic(n,d,L,sigma,k,alpha,delta,smearing):
    """
    Return transmission T based on analytical expressions for W and S
    """
    invLbAGS=(4.*sigma**2 / d**6 *(n*pi)**4 /
              k**2 * W_analytic(2*k,alpha,delta,smearing))
    invLbSGS=(0.5*(sigma/d*pi*n)**4 / k**2 *
              (1./3. + 1./(pi*n)**2)**2 * S_analytic(2*k,alpha,delta,smearing))
    return np.exp(-L*(invLbAGS+invLbSGS)), invLbAGS, invLbSGS

def T_numeric(n,d,L,sigma,kOld,kNew,W,S):
    """
    Return transmission T based on numerical expressions for W and S
    """
    # W and S are symmetric w.r.t. k - use only positive k-values for interpolation:
    kOld, W, S = [ x[len(x)/2:] for x in (kOld, W, S)]
    # Find interpolating values at k-values kNew:
    W, S = [ interp1d(kOld,x, kind="cubic")(2.*kNew) for x in (W, S)]
    invLbAGS=4.*sigma**2 / d**6 *(n*pi)**4 / kNew**2 * W
    invLbSGS=0.5*(sigma/d*pi*n)**4 / kNew**2 * (1./3. + 1./(pi*n)**2)**2 * S
    return np.exp(-L*(invLbAGS+invLbSGS)), invLbAGS, invLbSGS
