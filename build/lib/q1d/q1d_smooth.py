#!/usr/bin/python
# -*- coding: utf-8 -*-

# Sonst ist a/b=0 für a=2 und b=1
# aber a/b=0.5 für a=2. und b=1
from __future__ import division
import numpy as np
import q1d_utils



def rfftfreq(n_rfft,L):
    """ np.fft.fftfreq pendant for rfft: only generate positive  frequencies 
    Example:
    rfft  = np.fft.rfft(data**2)*Deltax
    rfftfreq= np.pi*q1d.rfftfreq(rfft.size,L) """
    n = n_rfft*2
    return np.fft.fftfreq(n,0.5*L/n)[0:n_rfft]

def symarray(c):
    """ Shift c[0] to the center of the array and mirror it """
    return np.concatenate((c[:0:-1],c))

def symaxis(c):
    """ Shift c[0] to the center of the array and mirror it, multiplied by -1"""
    return np.concatenate((-1*c[:0:-1],c))

# Aktuell beste Funktion für T:

def TasFTddCSquare(data,L,uptok=400):
    c = q1d_utils.NormAutoCorr(data)
    x,ddc = q1d_utils.Diff2(c,-L,L)
    # q1d_utils.Diff2 liefert array von -L bis L, daher *2
    Deltax=2*L/ddc.size
    kachse = np.linspace(0,uptok,uptok)
    i = 1j
    T = [ (np.exp(-i*k*x)*ddc**2).sum()*Deltax for k in kachse ]
    ImagReRatio = max(np.imag(T)/np.real(T))
    if ImagReRatio < 1e-5:
        T=np.real(T)
    else:
        print('TasFTddCSquare liefert imaginäres Ergebnis! Ratio Imag/Re ='+str(ImagReRatio))
    return [ kachse,T ]

def Tk(data,L,uptok=400):
    return TasFTddCSquare(data,L,uptok)

# Aktuell beste Funktion für W:

def WasFTC(data,L,uptok=400,Deltak=1):
    c = q1d_utils.NormAutoCorr(data)
    x = np.linspace(-L,L,c.size)
    Deltax=np.diff(x).mean() # = 2*L/c.size sonst
    assert 2*L/c.size - np.diff(x).mean() < 1e-6, 'Delta x falsch'
    kachse = np.linspace(0,uptok,int(uptok/Deltak))
    i = 1j
    W = [ (np.exp(-i*k*x)*c).sum()*Deltax for k in kachse ]
    ImagReRatio = max(np.imag(W)/np.real(W))
    if ImagReRatio > 1e-5:
        print('WasFTC liefert imaginäres Ergebnis!i') 
        print('Wir verwerfen das! Ratio war Im/Re ='+str(ImagReRatio))
    W=np.real(W)
    return [ kachse,W ]

def Wk(data,L,uptok=400):
    return WasFTC(data,L,uptok)

def WandT2k(data,L,upto=400):
    """ Calculate both W and T and scale to T(2k) and W(2k) """ 
    if L<0.1:
   	print 'Achtung L>0.1! L und sigma möglicherweise vertauscht?'
    kachseT,T = Tk(data,L,upto)
    kachseW,W = Wk(data,L,upto)
    assert (kachseW == kachseT).all(), 'Achtung Achsen falsch!'
    # Drop every other element
    W = W[::2]
    T = T[::2]
    assert (W.size == T.size), 'Arraylänge falsch'
    kachse = kachseW[:W.size]
    return kachse,W,T

def srbAloc(n,m):
    ret = 2.* np.pi**4 * n**2 * m**2 * np.cos(np.pi * (n - m)/2.)**4
    if ret < 1e-3:
        ret = 0
    return ret

def srbBloc(n,m):
    if n==m:
        ret = 0.25 * ((np.pi**4 * n**4)/9. + 2./3. * (np.pi * n)**2 + 1)
    else:
        ret = 64 * (n**2 * m**2 * (n**2 + m**2)**2)/(n**2 - m**2)**4 * np.cos((np.pi * (n  - m)/2)**4)
    if ret < 1e-3:
        ret = 0
    return ret

def srbAb(n,m):
    return 2.*srbAloc(n,m)

def srbBb(n,m):
    return 2.*srbBloc(n,m)

def GeneralOneOverLb(n,kachse,W,T,sigma,d=0.1):
    moden = np.arange(7)# TODO:: Prüfe wieviele es tatsächlich gibt...
    # m,n=0 existiert interessiert mich aber nicht!
    MsrbA = [ srbAb(n,m) for m in moden ]
    MsrbA[0] = 0

    MsrbB = [ srbBb(n,m) for m in moden ]
    MsrbB[0] = 0

    Mkachse = [ kn2km(kachse,n,m) for m in moden ]
    Mkachse[0] = Mkachse[0]*0

    Mk = np.array([ [ abs(kachse - k).argmin() for k in Mkachse[m] ] for m in moden ])
    Mtfwd = Mk[:]
    Mtbwd = Mk[:]
    # 
    # kn+km ist backward (alle Moden tragen bei)
    Mtbwd[n,m]= sigma**2/d**6 * srbA(n,m) * np.interp(Mk[n]+Mk[m], kachse, W)/kachse**2 \
              + sigma**4/d**4 * srbB(n,m) * np.interp(Mk[n]+Mk[m], kachse, T)/kachse**2 

    # kn-km ist forward  (nur n!=m trägt bei)
    Mtfwd[n,m]= sigma**2/d**6 * srbA(n,m) * np.interp(abs(Mk[n]-Mk[m]), kachse, W)/kachse**2 \
              + sigma**4/d**4 * srbB(n,m) * np.interp(abs(Mk[n]-Mk[m]), kachse, T)/kachse**2 
    # daher müssen die offdiagonalen terme gelöscht werden!
    Mtfwd = Mtfwd - Mtfwd*np.diag(np.ones(moden))

    return Widx
    #return sigma**2/d**6 * srbA(n,n) * W/kachse**2 + sigma**4/d**4 * srbB(n,n) 
    # TODO: W muss noch entsprechend verschoben werden

def OneOverLb(n,kachse,W,T,sigma,d=0.1):
    moden = [1,2,3]
    Mkachse = [ kn2km(kachse,n,m) for m in moden ]
    Mkachse[0] = Mkachse[0]*0

    Mk = np.array([ [ abs(kachse - k).argmin() for k in Mkachse[m] ] for m in moden ])
    Mtfwd = Mk[:]
    Mtbwd = Mk[:]
    # 
    # kn+km ist backward (alle Moden tragen bei)
    Mtbwd[n,m]= sigma**2/d**6 * srbA(n,m) * np.interp(Mk[n]+Mk[m], kachse, W)/kachse**2 \
              + sigma**4/d**4 * srbB(n,m) * np.interp(Mk[n]+Mk[m], kachse, T)/kachse**2 

    # kn-km ist forward  (nur n!=m trägt bei)
    Mtfwd[n,m]= sigma**2/d**6 * srbA(n,m) * np.interp(abs(Mk[n]-Mk[m]), kachse, W)/kachse**2 \
              + sigma**4/d**4 * srbB(n,m) * np.interp(abs(Mk[n]-Mk[m]), kachse, T)/kachse**2 
    # daher müssen die offdiagonalen terme gelöscht werden!
    Mtfwd = Mtfwd - Mtfwd*np.diag(np.ones(moden))

    return Widx
    #return sigma**2/d**6 * srbA(n,n) * W/kachse**2 + sigma**4/d**4 * srbB(n,n) 
    # TODO: W muss noch entsprechend verschoben werden



def t(n,kachse,W,T,sigma,L,d=0.1):
   return np.exp(-L*OneOverLb(n,kachse,W,T,sigma,d))

# Lbsm = L_b, single mode
def OneOverLbsm(n,kachse,W,T,sigma,d=0.1):
    return sigma**2/d**6 * srbAb(n,n) * W/kachse**2 + sigma**4/d**4 * srbBb(n,n) * T/kachse**2 

def tsm(n,kachse,W,T,sigma,L,d=0.1):
   return np.exp(-L*OneOverLbsm(n,kachse,W,T,sigma,d))
   # ist gleich exp(-L/Lb)=exp(-2L/Lloc)
