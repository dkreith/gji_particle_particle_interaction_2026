import numpy as np
import matplotlib.pyplot as plt
from functions import load_spec, scaleMaxwellGarnett

def calc_nu(a, d):
    Vk = (a**2*np.pi)/2
    Vd = (2*a+d)**2/2
    nu = Vk/Vd
    return nu

plt.rcParams.update({'font.size': 14})

filenames = ["1e-8", "2e-8", "5e-8", "1e-7", "2e-7", "5e-7", "1e-6", "2e-6",
             "5e-6", "1e-5", "2e-5"]
legend = ["$d = 1 $\u00D7$ 10^{-8}$ m", "$d = 2 $\u00D7$ 10^{-8}$ m",
          "$d = 5 $\u00D7$ 10^{-8}$ m", "$d = 1 $\u00D7$ 10^{-7}$ m",
          "$d = 2 $\u00D7$ 10^{-7}$ m", "$d = 5 $\u00D7$ 10^{-7}$ m",
          "$d = 1 $\u00D7$ 10^{-6}$ m", "$d = 2 $\u00D7$ 10^{-6}$ m",
          "$d = 5 $\u00D7$ 10^{-6}$ m", "$d = 1 $\u00D7$ 10^{-5}$ m",
          "$d = 2 $\u00D7$ 10^{-5}$ m"]
legend.append("rescaled")

dat = []

for filename in filenames:
    lspec = load_spec("Data/infinite_cluster", "d_"+filename+".txt")
    w = 2*np.pi*lspec[0][0:-9]
    sig = lspec[1][0:-9]
    nu = calc_nu(5e-6, float(filename))
    sig_rs = scaleMaxwellGarnett(nu, sig, 0.1, dim="2d")
    dat.append([w, sig_rs])
    
lspec0 = load_spec("Data/one_sphere/2d", "p_50.txt")
lspec0_load = scaleMaxwellGarnett(nu, lspec0[1], 0.1, dim="2d")
dat.append([lspec0[0], lspec0_load])
  
c0 = 1
F = 96485.33
mu = 5e-8
sig0 = 2*mu*F*c0
epsr = 80
eps0 = 8.85e-12  

col = plt.cm.viridis(np.linspace(0.1, 1, 11))
ii = 0

fig, ax = plt.subplots(figsize=(9,6))

for spec in dat:
    if ii == 11:
        ax.semilogx(spec[0], 1e4*spec[1].imag, "k--", linewidth=2)
    else:
        ax.semilogx(spec[0], 1e4*spec[1].imag, color=col[ii], linewidth=2)
    ii += 1

ax.legend(legend, loc="upper left")
ax.set_xlabel("$\u03C9$ [rad/s]")
ax.set_ylabel("$\u03C3''/\u03C3_0$ [-]")
ax.set_title("1e-4", loc="left", fontsize=14)
fig.tight_layout()
fig.savefig("Figures/infinite_cluster.pdf", dpi=300, bbox_inches="tight")