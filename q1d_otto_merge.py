#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import division
import pdb
import numpy as np

# Sonst ist a/b=0 für a=2 und b=1
# aber a/b=0.5 für a=2. und b=1

c = 2.99e8 # m/s

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


def HoleDaten(NrRepetitions,randomize=False):
    L = 30*5*2./1000. # 30 Pins a 5mm, jeder Pin doppel (*3) in Meter (/1000)
    if randomize:
        reihe = 1
    else:
        reihe = 0
    data = np.loadtxt("/home/otto/dat/corred_and_shuffeled.dat")[:,reihe]
    data = (0.8 + data*0.2 )/100.; # in mm
    data = ZeroMean(data)
    data,sigma = NormVar(data)
    data = np.array([ [ a ] * NrRepetitions for a in data]).flatten()
    return [data,sigma,L]


def HoleCNCDaten():
       L=0.8
       file = '/home/otto/dat/CNC/px260.dat'
       achse,data=np.loadtxt(file,unpack=True)
       data = ZeroMean(data)/1000. # in Meter
       data,sigma = NormVar(data)
       return [data,sigma,L]

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
    achse = np.linspace(vonL,bisL,x.size)
    dx=(bisL-vonL)/x.size
    print 'Diff: Boundary discretized with dx='+str(dx*1000)+' mm'
    return [achse[0:-1],np.diff(x)/np.diff(achse)]


def NormAutoCorr(x):
    result = np.correlate(x, x, mode='full')
    result=result/result.max()
    return result

def WasPowerSpectrum(data,L):
    c   = NormAutoCorr(data)
    ic  = np.fft.fftshift(c)
    Deltax  = 2*L/c.size
    fft = np.fft.fft(ic)*Deltax
	 # Falsch, aber äquivalent zur neuen Version:
    # # Pi wegen Exp(i\pi)ax)->Exp(ikx)
    # fftachse=np.pi*np.fft.fftfreq(ic.size,L/ic.size)

    # 2 Pi wegen Exp(i*2*pi*a*x)->Exp(ikx)
    fftachse=2*np.pi*np.fft.fftfreq(ic.size,Deltax)
    return [fftachse,fft]

def WasFTC(data,L,uptok=400,Deltak=1):
    c = NormAutoCorr(data)
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
def Tk(data,L,uptok=400):
    return TasFTddCSquare(data,L,uptok)

# Aktuell beste Funktion für W:
def Wk(data,L,uptok=400):
    return WasFTC(data,L,uptok)


def TasFTddCSquare(data,L,uptok=400):
    c = NormAutoCorr(data)
    x,ddc = Diff2(c,-L,L)
    # Diff2 liefert array von -L bis L, daher *2
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

def TasRFTddCSquare(data,L):
    c = NormAutoCorr(data)
    x,ddc = Diff2(c,-L,L)
    #Das ganze als rfft
    # 2.0 * rfft = fft 
    # weil fft doppelt so lang ist.
    # bzw. rfft nur über das halbe Array geht
    data = ddc[ddc.size/2:]
    Deltax=2*L/ddc.size
    rfft  = 2.0*np.fft.rfft(data**2)*Deltax 
    # Pi wegen Exp(-i\pi)ax)->Exp(-ikx)
	 print "ACHTUNG: Hier muss noch was korrigiert werden!, *np.pi ist falsch
	 und L als letzter Parameter auch " # 2 Pi wegen Exp(i*2*pi*a*x)->Exp(ikx)
    freq= np.pi*rfftfreq(rfft.size,L)
    return freq,rfft

def TasFFTddCSquare(data,L):
    c = NormAutoCorr(data)
    x,ddc = Diff2(c,-L,L)
    iddc = np.fft.fftshift(ddc)
    Deltax=2*L/ddc.size
    fft   = np.fft.fft(iddc**2)*Deltax 
    # 2Pi wegen Exp(i2\pi*ax)->Exp(ikx)
    fftfreq=2*np.pi*np.fft.fftfreq(iddc.size,Deltax)
    return fftfreq,fft

# DIE FUNKTIONIEREN NOCH NICHT
def BROKENTasConv(sdata,L):
    xi_k=np.fft.rfft(sdata)
    k=rfftfreq(xi_k.size,L)# b.size=sdata.size+1
    c=abs(k*xi_k)**2
    c=symarray(c)
    k=symaxis(k)
    return k,np.convolve(c,c,mode='same')

def BROKENTasFTC1Square(x,L,window):
    achse,diffx = SmoothDiff(x,L,window)
    cdiffx = np.correlate(diffx, diffx, mode='full')
    # N1 = (diffx**2).sum() # Sum over \xi'(x)^2 == np.correlate(diffx, diffx, mode='full')
    # N1 = (diffx**2).sum()*L/diffx.size # Integral over \xi'(x)^2 
    #    ==  L/diffx.size*np.correlate(diffx, diffx, mode='full')
    dx = L/diffx.size
    cdiffx =  cdiffx*dx
    # np.correlate reutrns centered array, i.e., f(0) is at N/2
    # fft needs f(0) to be at 0, and positive freq components at 1..N/2
    # negative freq components at N/2..N. fftshift fixes this
    cdiffx   = np.fft.fftshift(cdiffx)
    fftachse = np.fft.fftfreq(cdiffx.size,dx)
    fft = np.fft.fft(cdiffx**2)/L**2 # L=N0
    return [fftachse,fft]

def WandT2k(data,L,upto=400):
    assert L>0.1, 'L und sigma möglicherweise vertauscht'
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
   # ist gleich exp(-L/Lb)=exp(-1L/Lloc)
