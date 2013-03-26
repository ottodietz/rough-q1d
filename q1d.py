#!/usr/bin/python
# -*- coding: utf-8 -*-

# Sonst ist a/b=0 für a=2 und b=1
# aber a/b=0.5 für a=2. und b=1
from __future__ import division

import numpy as np

c = 2.99792458e8 # m/s

def frq2k(nu):
    return 2. * np.pi *  nu /c;#(*1/s -> 1/m*)

def k2frq(k):
    return (c*k)/(2.*np.pi) #(*1/m -> 1/s*)

def frq2wvl(nu):
    return c/nu

def wvl2frq(wvl):
    return c/wvl#(* nu=c/lambda*)

def kx2frq(kx,  mode, d = 0.1) :
    return np.sqrt(k2frq(kx)**2 + (c*mode/(2.*d))**2)

def frq2kx(frq, mode, d = 0.1) :
    return frq2k(np.sqrt(frq**2 - (c*mode/(2.*d))**2))

def kn2km(kachse,kfrom,kto):
    return frq2kx(kx2frq(kachse,kfrom),kto)

def modeopensatfrq(n,d=0.1) :
    return wvl2frq(2.*d/n)

def NormVar(x):
    sigma = np.sqrt(x.var())
    x = x / sigma 
    return [x,sigma]

def ZeroMean(x):
    x = x - x.mean()
    return x


def Smooth(data,WINDOW): 
    assert WINDOW > 1, 'window size <= 1, yields empty result'
    weightings = np.ones(WINDOW) / WINDOW
    # Kann ich hier statt dessen einfach 'same' benutzen?
    return np.convolve(data, weightings)[WINDOW-1:-(WINDOW-1)]

def Diff2(x,vonL,bisL):
    """ return the second derivative of x """ 
    achse = np.linspace(vonL,bisL,x.size)[1:]
    #ddx = (np.diff(x)- 2*x[0:-1] + np.diff(x[::-1]))/(np.diff(achse)**2)
    fx=x[1:-1]#f(x)        =   - 2 3 4 -
    fm=x[:-2]#f(x-Deltax)  =     1 2 3 - -
    fp=x[2:]#f(x+Deltax)   = - - 3 4 5
    ddx = (fp+fm-2*fx)/(np.diff(achse)**2)
    dx=(bisL-vonL)/x.size
    print 'Diff2: Boundary discretized with dx='+str(dx*1000)+\
          ' mm. Zum Vergleich, np.mean(np.diff(achse))='+str(np.mean(np.diff(achse))*1000)+'mm'
    return achse[:-1], ddx

def Diff(x,vonL,bisL):
    """ return the derivative of x """ 
    achse = np.linspace(vonL,bisL,x.size)
    dx=(bisL-vonL)/x.size
    print 'Diff: Boundary discretized with dx='+str(dx*1000)+' mm'
    return [achse[0:-1],np.diff(x)/np.diff(achse)]


def NormAutoCorr(x):
	 """ correlate x with itself and normalize to C(0)=1 """
    result = np.correlate(x, x, mode='full')
    result=result/result.max()
    return result

