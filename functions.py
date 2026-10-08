import numpy as np
import scipy.interpolate as interp
import scipy.signal as sgnl

def load_spec(path, fname):
    dat = np.loadtxt(path+"/"+fname, skiprows=5, dtype=str)

    def make_comp(string):
        real = ""
        imag = ""
        isreal = True
        realexp = False
    
        for c in string:
            if isreal == True:
                if (c == "+" or (c == "-" and real != "")) and realexp==False:
                    imag += c
                    isreal = False
                else:
                    if c == "E":
                        realexp = True
                    real += c
                    if realexp == True and (c == "+" or c == "-"):
                        realexp = False
            elif isreal == False and c != "i":
                imag += c
        comp = float(real) + 1j * float(imag)
        return comp

    w = []
    sig = []

    for datpoint in dat:
        w.append(float(datpoint[0]))
        sig.append(make_comp(datpoint[1]))
    
    w = np.array(w)
    sig = np.array(sig)
    
    return w, sig

def debyeLength(epsr, T, C0, mod="concentration"):
    
    eps0 = 8.85e-12
    kB = 1.38e-23
    e = 1.602e-19
    
    if mod == "concentration":
        NA = 6.022e23
    elif mod == "density":
        NA = 1
    
    c0 = C0*NA
    
    lamb = np.sqrt(eps0*epsr*kB*T/2/c0/e**2)
    return lamb

def maxwellGarnett(sigi, sigm, nu, f=None, mod="norm", dim="3d"):
    
    if np.any(f == None):
        if dim == "3d":
            fmg = (sigi-sigm)/(2*sigm+sigi)
        elif dim == "2d":
            fmg = (sigi-sigm)/(sigm+sigi)
    else:
        fmg = f
    
    if dim == "3d":
        signorm = (1+2*nu*fmg)/(1-nu*fmg)
    elif dim == "2d":
        signorm = (1+nu*fmg)/(1-nu*fmg)
    sigmg = sigm*signorm
    
    if mod == "sig":
        return sigmg
    elif mod == "norm":
        return signorm
    elif mod == "f":
        return fmg
    
def scaleMaxwellGarnett(nu, sig, nu0, dim="3d"):
    if dim == "2d":
        f = 1/nu*(sig-1)/(sig+1)
    elif dim == "3d":
        f = 1/nu*(sig-1)/(sig+2)
    sig_sc = maxwellGarnett(0, 0, nu0, f=f, dim=dim)
    return sig_sc

def lyklema(sigi, sigm, nu, a, w, muS, SigS, C0, epsr, T, zeta,
            cmod="concentration", smod="norm", kmod=False):
    
    kB = 1.38e-23
    e = 1.6e-19
    
    DS = muS*kB*T/e
    
    if cmod == "concentration":
        NA = 6.022e23
    elif cmod == "density":
        NA = 1
        
    c0 = C0*NA
    
    lambdaD = debyeLength(epsr, T, C0, mod=cmod)
    M = 1+SigS/(2*e*c0*lambdaD*np.cosh((e*zeta)/(2*kB*T)))
    
    tauL = a**2/(2*DS*M)
    kappaL = (1j*w*tauL)/(1+1j*w*tauL)*muS*SigS
    f = (sigi+2*kappaL/a-sigm)/(2*sigm+sigi+2*kappaL/a)
    
    sig = maxwellGarnett(sigi, sigm, nu, f=f, mod=smod)
    
    if kmod:
        return kappaL
    else:
        return sig

def findmax(w, sigma, high=0):
    
    # Interpolation
    w_int = np.logspace(np.log10(np.amin(w)),
                        np.log10(np.amax(w)),1000)
    sig_int = interp.interp1d(w, sigma, kind='cubic', fill_value='extrapolate')
    
    sig_int_arr = sig_int(w_int)
    
    peaks, _ = sgnl.find_peaks(sig_int_arr)
    if len(peaks)==0:
        sig_max = np.nan
        tau = np.nan
    else:
        #if high == 0:
        if (sig_int_arr[peaks[0]] > sig_int_arr[peaks[len(peaks)-1]] and
            len(peaks) > 1):
            index = peaks[0]
        #elif high == 1:
        else:
            index = peaks[len(peaks)-1]
        w_max = w_int[index]
        tau = 1/w_max
    
        sig_max = sig_int_arr[index]
    
    return [sig_max, tau]