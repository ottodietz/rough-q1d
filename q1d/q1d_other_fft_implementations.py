
# Other implementations

def WasPowerSpectrum(data,L):
    c   = q1d_utils.NormAutoCorr(data)
    ic  = np.fft.fftshift(c)
    Deltax  = 2*L/c.size
    fft = np.fft.fft(ic)*Deltax
	 # Falsch, aber äquivalent zur neuen Version:
    # # Pi wegen Exp(i\pi)ax)->Exp(ikx)
    # fftachse=np.pi*np.fft.fftfreq(ic.size,L/ic.size)

    # 2 Pi wegen Exp(i*2*pi*a*x)->Exp(ikx)
    fftachse=2*np.pi*np.fft.fftfreq(ic.size,Deltax)
    return [fftachse,fft]


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
    print "ACHTUNG: Hier muss noch was korrigiert werden!, *np.pi ist falsch\
	 			und L als letzter Parameter auch " 
    # 2 Pi wegen Exp(i*2*pi*a*x)->Exp(ikx)
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
