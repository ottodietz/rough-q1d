#!/usr/bin/env python

import sys
import argparse
from math import *
import numpy as np
from numpy.fft import fft, fftshift
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

# ---------------------------------------------
# Parsing command-line arguments
# ---------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument("-m", "--mode", default="1",
                    help="set wire mode" )
parser.add_argument("-s", "--smearing", default="0.001", nargs="+", type=float, 
                    help="set smearing parameter (allows multiple arguments)" )
parser.add_argument("-a", "--alpha-data", default="",
                    help="set directory to pin data" )
parser.add_argument("-e", "--experimental-data", default="../Exp_data/",
                    help="set directory to experimental data" )
parser.add_argument("-n", "--numerical-data", default="../Num_data/",
                    help="set directory to numerical data" )
args = parser.parse_args()

mode = args.mode
smearing = args.smearing
if type(smearing) is not list:
    smearing = [ args.smearing ]
alpha_path = args.alpha_data
exp_path = args.experimental_data
num_path = args.numerical_data

# -----------------------------------------------------------------------------
# Constants
# -----------------------------------------------------------------------------
alphaFile = alpha_path + "corred_and_shuffeled.dat"  # pinheight data
expFiles = [ exp_path + x + "_mode" + mode + "_kx.dat" for x in "rscb", "rsrb" ]
numFiles = [ num_path + "T" + 2*mode + "_num_" + x + ".dat" for x in "c", "s" ]

pinSize = 0.005                             # pin size
delta = 10*pinSize                          # module width
wireLength = 50*delta                            # wire length L
wireWidth = 0.1                             # wire width d

samplingLength = (wireLength+400.*delta)    # sampling length
dx = .0001                                 # step-length in x-space
dk = 2.*pi/samplingLength                   # step-length in k-space
xRange = np.arange(-samplingLength/2,
                   samplingLength/2,dx)     # discretization x-axis
kRange = xRange*dk/dx                      # discretization k-axis
kRange = np.arange(0.01,600,0.01)
kRangeNew = np.arange(0.01,600,0.01)        # k-values for plot

# -----------------------------------------------------------------------------
# Read data files and generate pin heights
# -----------------------------------------------------------------------------
alpha = [ (x - x.mean()) for x in
          ( (0.8 + 0.2*y)/100 for y in np.loadtxt(alphaFile,unpack=True) ) ]
exp = [ np.loadtxt(i,unpack=True) for i in expFiles ]
num = [ np.loadtxt(i,unpack=True) for i in numFiles ]
# randomly chosen steps:
alpha_shuffle = np.random.uniform(-sqrt(3.)/100,sqrt(3.)/100, size=50)
alpha_corr = np.random.uniform(-sqrt(3.)/100,sqrt(3.)/100, size=50)
alpha = [ alpha_shuffle, alpha_corr ]

## -----------------------------------------------------------------------------
## Function definitions
## -----------------------------------------------------------------------------

### ALL FUNCTIONS MOVED TO q1d_step
# to make this script working:
from q1d_step import *

# -----------------------------------------------------------------------------
# Plot data: xi, W, S (analytic)
# -----------------------------------------------------------------------------
if False:
    # -------------------------------------------------------------------------
    # Plot data: T (theory, simulation, experiment)
    # -------------------------------------------------------------------------
    for alpha, exp, num, subplotindex, label in zip(alpha, exp, num, [1,2],
                                                    ["Correlated","Shuffled"]):
        plt.subplot(2,1,subplotindex)
        plt.semilogy(num[0],num[int(mode)+4],"r-o",
                     label=r"$T_{%s}$ simulation" % (2*mode))  
        plt.semilogy(exp[0],exp[1]*1000,"k-",
                     label=r"$T_{%s}$ experiment ($\times 10^3$)" % (2*mode))
        for s in smearing:
            T,_,_ = T_analytic(int(mode),wireWidth,wireLength,
                               1.,kRangeNew,alpha,delta,s)
            plt.semilogy(kRangeNew,T,
                         label=r"$T_{%s}$ ($\sigma$ = %g)" % (2*mode,s)) 
        plt.title(label)
        plt.xlim(0,400)
        plt.ylim(1e-7,1e1)
        plt.legend(loc=4)
        plt.ylabel(r"$T_{%s}$" % (2*mode),fontsize="x-large")
    plt.xlabel(r"$k_x$",fontsize="x-large")
    plt.show()

# -----------------------------------------------------------------------------
# Plot data: xi, W, S (numeric)
# -----------------------------------------------------------------------------
if True:
    # ---------------------------------------------
    # Fourier transform for roughness-height power spectrum W(k)
    # ---------------------------------------------
    xi = xi_corr, xi_shuffle = [ build_wire(xRange,a,delta,smearing) for a in alpha ]  # discretized wire profile xi(x)
    W_numeric_corr, W_numeric_shuffle = [ fourier(i) / wireLength for i in xi ]  # Fourier transform wire profile xi and shift data (centered zero frequency)
    W_analytic_corr, W_analytic_shuffle = [ W_analytic(kRange,a,delta,smearing[0]) for a in alpha  ]  # analytical W(k)
    
    # ---------------------------------------------
    # Fourier transform for square-gradient power spectrum S(k)
    # ---------------------------------------------
    xipsq = xipsq_corr, xipsq_shuffle = [ np.gradient(x/dx)**2 for x in xi ]  # discretized square gradient xi'(x)^2
    S_numeric_corr, S_numeric_shuffle = [ fourier(i) / wireLength / 2 for i in xipsq ]  # factor 1/2 from <V(x)V(x')> correlator definition
    S_analytic_corr, S_analytic_shuffle = [ S_analytic(kRange,a,delta,smearing[0]) for a in alpha  ]  # analytical S(k)
    
    # ---------------------------------------------
    # Transmission T(k), 1/L^(b,AGS), 1/L^(b,SGS)
    # ---------------------------------------------
    #T, invLbAGS, invLbSGS = T_numeric(3.,wireWidth,wireLength,1.,kRange,kRangeNew,W_numeric_shuffle,S_numeric_shuffle)
    
    # ---------------------------------------------
    # Plot data: xi, W, S
    # ---------------------------------------------
    p1 = plt.figure(figsize=(16,7))
    # compare step-profiles
    plt.subplot(311)
    plt.xlim(-0.2,0.2)
    plt.ylim(-0.15,0.15)
    plt.plot(xRange,-xi_shuffle+wireWidth,"r-",label=r"$\xi(x)$ shuffled")
    plt.plot(xRange,xi_shuffle,"r-")
    plt.plot(xRange,-xi_corr+wireWidth-0.05,"g-",label=r"$\xi(x)$ correlated")
    plt.plot(xRange,xi_corr-0.05,"g-")
    plt.legend(loc=8,ncol=2)
    # compare W power spectra
    plt.subplot(312)
    plt.xlim(0,120)
    plt.ylim(1e-12,1e-5)
    plt.semilogy(kRange/2,W_analytic_shuffle,"r-",label="W(k) shuffled")
    plt.semilogy(kRange/2,W_analytic_corr,"g-",label="W(k) correlated")
    plt.legend(loc=8,ncol=2)
    # compare S power spectra
    plt.subplot(313)
    plt.xlabel("k",fontsize=14)
    plt.xlim(0,120)
    plt.ylim(1e-6,1e2)
    plt.semilogy(kRange/2,S_analytic_shuffle,"r-",label="S(k) shuffled")
    plt.semilogy(kRange/2,S_analytic_corr,"g-",label="S(k) correlated")
    plt.legend(loc=8,ncol=2)

    plt.show()
