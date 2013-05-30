import numpy as np

def srbAloc(n,m):
    """ Coefficient for L_{loc} """
    ret = 2.* np.pi**4 * n**2 * m**2 * np.cos(np.pi * (n - m)/2.)**4
    if ret < 1e-3:
        ret = 0
    return ret

def srbBloc(n,m):
    """ Coefficient for L_{loc} """
    if n==m:
        ret = 0.25 * ((np.pi**4 * n**4)/9. + 2./3. * (np.pi * n)**2 + 1)
    else:
        ret = 64 * (n**2 * m**2 * (n**2 + m**2)**2)/(n**2 - m**2)**4 * np.cos((np.pi * (n  - m)/2)**4)
    if ret < 1e-3:
        ret = 0
    return ret

def srbAb(n,m):
    """ Coefficient for L_{backscattering} """
    return 2.*srbAloc(n,m)

def srbBb(n,m):
    """ Coefficient for L_{backscattering} """
    return 2.*srbBloc(n,m)

# Lbsm: Lb = L_b (backscattering length), sm = single mode
def OneOverLbsm(n,kachse,W,T,sigma,d=0.1):
    """ 1/L_{backscattering} for single mode (sm) """
    return sigma**2/d**6 * srbAb(n,n) * W/kachse**2 + sigma**4/d**4 * srbBb(n,n) * T/kachse**2 

def tsm(n,kachse,W,T,sigma,L,d=0.1):
    """ t=exp(-L/L_b) for single mode (sm) """
    return np.exp(-L*OneOverLbsm(n,kachse,W,T,sigma,d)) 

###
### Under construction
###
def GeneralOneOverLb(n,kachse,W,T,sigma,d=0.1):
    moden = np.arange(7)# TODO: Wievile gibt es wirklich?
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

    # kn-km ist forward  (nur n!=m tragen bei)
    Mtfwd[n,m]= sigma**2/d**6 * srbA(n,m) * np.interp(abs(Mk[n]-Mk[m]), kachse, W)/kachse**2 \
              + sigma**4/d**4 * srbB(n,m) * np.interp(abs(Mk[n]-Mk[m]), kachse, T)/kachse**2 
    # daher muessen die offdiagonalen terme geloescht werden!
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

    # kn-km ist forward  (nur n!=m tragen bei)
    Mtfwd[n,m]= sigma**2/d**6 * srbA(n,m) * np.interp(abs(Mk[n]-Mk[m]), kachse, W)/kachse**2 \
              + sigma**4/d**4 * srbB(n,m) * np.interp(abs(Mk[n]-Mk[m]), kachse, T)/kachse**2 
    # daher muessen die offdiagonalen terme geloescht werden!
    Mtfwd = Mtfwd - Mtfwd*np.diag(np.ones(moden))

    return Widx
    #return sigma**2/d**6 * srbA(n,n) * W/kachse**2 + sigma**4/d**4 * srbB(n,n) 
    # TODO: W muss noch entsprechend verschoben werden

def t(n,kachse,W,T,sigma,L,d=0.1):
   return np.exp(-L*OneOverLb(n,kachse,W,T,sigma,d))





