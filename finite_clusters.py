import numpy as np
import matplotlib.pyplot as plt
from functions import load_spec, scaleMaxwellGarnett

Vs = np.pi*(5e-6)**2
Vd = (5e-5)**2
nu = Vs/Vd

plt.rcParams.update({'font.size': 15})

filenames = ["2", "3", "4", "5", "6", "7", "8", "9", "10"]

dat = []

for filename in filenames:
    lspec = load_spec("Data/finite_cluster",
                      "d_1e-7_p_0,5_size_"+filename+".txt")
    if filename == "4" or filename == "8":
        lspec_load = scaleMaxwellGarnett(0.101, lspec[1], 0.1, dim="2d")
    elif filename == "7":
        lspec_load = scaleMaxwellGarnett(0.102, lspec[1], 0.1, dim="2d")
    else:
        lspec_load = lspec[1]
    dat.append([lspec[0], lspec_load])

lspec0 = load_spec("Data/one_sphere/2d", "p_50.txt")
lspec0_load = scaleMaxwellGarnett(nu, lspec0[1], 0.1, dim="2d")
dat.append([lspec0[0], lspec0_load])

lspec2_0 = load_spec("Data/angle_dependence/d_1e-7/two_spheres_p_50",
                     "d_1e-7_phi_0.txt")
lspec2_0_load = scaleMaxwellGarnett(nu*2, lspec2_0[1], 0.1, dim="2d")
dat.append([lspec2_0[0], lspec2_0_load])

lspec2_90 = load_spec("Data/angle_dependence/d_1e-7/two_spheres_p_50",
                      "d_1e-7_phi_90.txt")
lspec2_90_load = scaleMaxwellGarnett(nu*2, lspec2_90[1], 0.1, dim="2d")
dat.append([lspec2_90[0], lspec2_90_load])

c0 = 1
F = 96485.33
mu = 5e-8
sig0 = 2*mu*F*c0
epsr = 80
eps0 = 8.85e-12

col = plt.cm.viridis(np.linspace(0.1, 1, 9))
ii = 0

fig, ax = plt.subplots(figsize=(9,8.5))
for spec in dat:
    if ii == 9:
        ax.semilogx(spec[0], 1e5*spec[1].imag, "--k", linewidth=2)
    elif ii == 10:
        ax.semilogx(spec[0], 1e5*spec[1].imag, "-.", color="gray",
                    linewidth=2)
    elif ii == 11:
        ax.semilogx(spec[0], 1e5*spec[1].imag, ":", color="slategray",
                    linewidth=2)
    else:
        ax.semilogx(2*np.pi*spec[0],1e5* spec[1].imag, linewidth=2, color=col[ii])
    ii += 1
    
leg = []
for size in filenames:
    leg.append(size+"$\u00D7$"+size+" cluster")
leg.append("rescaled")
leg.append("2 part. ($\u03B1 = 0$°)")
leg.append("2 part. ($\u03B1 = 90$°)")
ax.legend(leg, loc="upper left")
ax.set_xlabel("$\u03C9$ [rad/s]")
ax.set_ylabel("$\u03C3''/\u03C3_0$ [-]")
ax.set_title("1e-5", loc="left", fontsize=15)
fig.tight_layout()
fig.savefig("Figures/finite_cluster.pdf", dpi=300, bbox_inches="tight")