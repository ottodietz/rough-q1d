import sys
import getopt
import argparse
import numpy as np
from numpy.fft import fft, fftshift
from scipy.interpolate import interp1d
from math import *
import matplotlib.pyplot as plt

# ---------------------------------------------
# Parsing command-line arguments
# ---------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument("-m", "--mode", default="1",
                    help="set wire mode" )
parser.add_argument("-s", "--smearing", default="0.001", type=float, nargs='+',
                    help="set smearing parameter (allows multiple arguments)" )
parser.add_argument("-e", "--experimental-data", default="../Exp_data/",
                    help="set directory to experimental data" )
parser.add_argument("-n", "--numerical-data", default="../Num_data/",
                    help="set directory to numerical" )
args = parser.parse_args()

mode = args.mode
smearing = args.smearing
exp_path = args.experimental_data
num_path = args.numerical_data

# ---------------------------------------------
# Constants
# ---------------------------------------------
alphaFile="corred_and_shuffeled.dat"  # input file for boundary data
expFile_corr, expFile_shuffle = [ exp_path + x + "_mode" + mode + "_kx.dat" for x in "rscb", "rsrb" ]
numFile_corr, numFile_shuffle = [ num_path + "T" + 2*mode + "_num_" + x + ".dat" for x in "c", "s" ]
#
delta = 0.01                            # module width
wireLength = 0.3                        # wire length L
wireWidth = 0.1                         # wire width d
#
sampleLength = (wireLength+400.*delta)  # sampling length
dx = .00005                             # step-length in x-space
dk = 2.*pi/sampleLength                 # step-length in k-space
xRange = np.arange(-sampleLength/2.,sampleLength/2.,dx)   # discretization of x-axis
kRange = xRange*dk/dx                                     # discretization of k-axis
kRangeNew = np.arange(0.01,200,0.25)                      # k-values for plot: (0,200)

# ---------------------------------------------
# Read data files and generate pin heights
# ---------------------------------------------
alpha = alpha_corr, alpha_shuffle = [ (x - x.mean()) for x in ( (0.8 + 0.2*y)/100 for
                                            y in np.loadtxt(alphaFile,unpack=True) ) ]
exp = exp_corr, exp_shuffle = [ np.loadtxt(i,unpack=True) for
                                i in (expFile_corr,expFile_shuffle)]
num = num_corr, num_shuffle = [ np.loadtxt(i,unpack=True) for
                                i in (numFile_corr,numFile_shuffle)]
#alpha_shuffle = np.random.uniform(-sqrt(3.),sqrt(3.), size=alpha_shuffle.size) # randomly chosen steps

# ---------------------------------------------
# Function definitions
# ---------------------------------------------
def build_wire(x,alpha,delta,smearing):
    "Return sum_(n=-N,N) alpha_n step(x-n*Delta)"
    
    def step(x,delta,smearing):
        "Return smoothed step-function Fermi(x-Delta)-Fermi(x)"
        return 1./(1.+np.exp((x-delta)/smearing))-1./(1.+np.exp(x/smearing))
    
    module_number = alpha.size
    temp_wire = np.zeros(x.size)
    for n in np.arange(module_number):
        temp_wire = ( temp_wire + alpha[n] *
                        step(x - (n-module_number/2.)*delta,delta,smearing) )
    return temp_wire

def fourier(data):
    "Return Fourier transform of input data"
    return fftshift(abs(fft(data)*dx))**2

def W_analytic(k,alpha,delta,smearing):
    "Return roughness-height power spectrum W(k)"
    module_index = np.arange(alpha.size)
    # Take correlation of alphas into account:
    omega = ([ np.exp(-1j*n*k*delta)*alpha[n] for n in module_index] )
    omega = np.sum(omega,axis=0)
    omega = abs(omega)**2 / alpha.size
    return ( 1./delta * (2.*pi*smearing*
                np.sinh(k*pi*smearing)**(-1)*np.sin(k*delta/2.))**2 ) * omega

def S_analytic(k,alpha,delta,smearing):
    "Return roughness-height power spectrum S(k)"
    module_index = np.arange(alpha.size)
    a = np.concatenate([ [0],alpha,[0] ])  # a[n-1] and a[n+1] = 0 for n=N and n=0
    # Take correlation of alphas into account:
    omega = ([ np.exp(-1j*n*k*delta)*(a[n]*(a[n]-a[n+1])*
                np.exp(-1j*k*delta) + a[n]*(a[n]-a[n-1])) for n in module_index] )
    omega = np.sum(omega,axis=0)  
    omega = abs(omega)**2 / alpha.size 
    return ( 1./delta / 72. * (k*pi*(1.+k**2*smearing**2)*
                            np.sinh(k*pi*smearing)**(-1))**2 ) * omega

def T_analytic(n,d,L,sigma,k,alpha,delta,smearing):
    "Return transmission T based on analytical expressions for W and S"
    invLbAGS=( 4.*sigma**2 / d**6 *(n*pi)**4 /
                k**2 * W_analytic(2*k,alpha,delta,smearing) )
    invLbSGS=( 0.5*(sigma/d*pi*n)**4 / k**2 *
                (1./3. + 1./(pi*n)**2)**2 * S_analytic(2*k,alpha,delta,smearing) )
    return np.exp(-L*(invLbAGS+invLbSGS)), invLbAGS, invLbSGS

def T_numeric(n,d,L,sigma,kOld,kNew,W,S):
    "Return transmission T based on numerical expressions for W and S"
    # W and S are symmetric w.r.t. k - use only positive k-values for interpolation:
    kOld, W, S = [ x[len(x)/2:] for x in (kOld, W, S)]
    # Find interpolating values at k-values kNew:
    W, S = [ interp1d(kOld,x, kind="cubic")(2.*kNew) for x in (W, S) ]
    invLbAGS=4.*sigma**2 / d**6 *(n*pi)**4 / kNew**2 * W
    invLbSGS=0.5*(sigma/d*pi*n)**4 / kNew**2 * (1./3. + 1./(pi*n)**2)**2 * S
    return np.exp(-L*(invLbAGS+invLbSGS)), invLbAGS, invLbSGS

# ---------------------------------------------
# Plot data: xi, W, S (analytic)
# ---------------------------------------------
if True:
    # ---------------------------------------------
    # Plot data: T (theory, simulation, experiment)
    # ---------------------------------------------
    for alpha, exp, num, subplotindex, label in zip( alpha, exp, num, [1,2], ["Correlated","Shuffled"] ):
        plt.subplot(2,1,subplotindex)
        plt.semilogy(num[0],num[int(mode)+4],"r-o",label=r"$T_{%s}$ simulation" % (2*mode) )  # simulation results
        plt.semilogy(exp[0],exp[1]*1000,"k-",label=r"$T_{%s}$ experiment ($\times 10^3$)" % (2*mode) )  # experimental results
        for s in smearing:
            T,_,_ = T_analytic(int(mode),wireWidth,wireLength,1.,kRangeNew,alpha,delta,s)
            plt.semilogy(kRangeNew,T,label=r"$T_{%s}$ ($\sigma$ = %g)" % (2*mode,s) )  # theory
        plt.title(label)
        plt.xlim(0,100)
        plt.ylim(1e-7,1e1)
        plt.legend(loc=4)
        plt.ylabel(r"$T_{%s}$" % (2*mode),fontsize="x-large")
    plt.xlabel(r"$k_x$",fontsize="x-large")
    plt.show()

# ---------------------------------------------
# Plot data: xi, W, S (numeric)
# ---------------------------------------------
if False:
    # ---------------------------------------------
    # Fourier transform for roughness-height power spectrum W(k)
    # ---------------------------------------------
    xi = xi_corr, xi_shuffle = [ build_wire(xRange,a,delta,smearing) for a in alpha ]  # discretized wire profile xi(x)
    W_numeric_corr, W_numeric_shuffle = [ fourier(i) / wireLength for i in xi ]  # Fourier transform wire profile xi and shift data (centered zero frequency)
    W_analytic_corr, W_analytic_shuffle = [ W_analytic(kRange,a,delta,smearing) for a in alpha  ]  # analytical W(k)
    
    # ---------------------------------------------
    # Fourier transform for square-gradient power spectrum S(k)
    # ---------------------------------------------
    xipsq = xipsq_corr, xipsq_shuffle = [ np.gradient(x/dx)**2 for x in xi ]  # discretized square gradient xi'(x)^2
    S_numeric_corr, S_numeric_shuffle = [ fourier(i) / wireLength / 2 for i in xipsq ]  # factor 1/2 from <V(x)V(x')> correlator definition
    S_analytic_corr, S_analytic_shuffle = [ S_analytic(kRange,a,delta,smearing) for a in alpha  ]  # analytical S(k)
    
    # ---------------------------------------------
    # Transmission T(k), 1/L^(b,AGS), 1/L^(b,SGS)
    # ---------------------------------------------
    #T, invLbAGS, invLbSGS = T_numeric(3.,wireWidth,wireLength,1.,kRange,kRangeNew,W_numeric_shuffle,S_numeric_shuffle)
    
    # ---------------------------------------------
    # Plot data: xi, W, S
    # ---------------------------------------------
    p1 = plt.figure(figsize=(16,10))
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
    plt.xlim(0,300)
    plt.ylim(1e-12,1e-5)
    plt.semilogy(kRange,W_analytic_shuffle,"r-",label="W(k) shuffled")
    plt.semilogy(kRange,W_analytic_corr,"g-",label="W(k) correlated")
    plt.legend(loc=8,ncol=2)
    # compare S power spectra
    plt.subplot(313)
    plt.xlabel("k",fontsize=14)
    plt.xlim(0,300)
    plt.ylim(1e-6,1e2)
    plt.semilogy(kRange,S_analytic_shuffle,"r-",label="S(k) shuffled")
    plt.semilogy(kRange,S_analytic_corr,"g-",label="S(k) correlated")
    plt.legend(loc=8,ncol=2)

    plt.show()