import numpy as np
import matplotlib.pyplot as plt
from functions import load_spec, findmax, scaleMaxwellGarnett

Vs_3 = 4/3*np.pi*(5e-6)**3
Vd_3 = np.pi*(5e-5)**2*(2*5e-5)
nu_3 = Vs_3/Vd_3
nu0_3 = 2*nu_3

Vs_2 = np.pi*(5e-6)**2
Vd_2 = (2*5e-5)**2
nu_2 = Vs_2/Vd_2
nu0_2 = 2*nu_2

filenames = ["1e-8", "15e-9", "2e-8", "3e-8", "5e-8", "7e-8",
             "1e-7", "15e-8", "2e-7", "3e-7", "5e-7", "7e-7",
             "1e-6", "15e-7", "2e-6", "3e-6", "5e-6", "7e-6",
             "1e-5", "15e-6", "2e-5", "3e-5", "5e-5", "7e-5",
             "1e-4", "15e-5", "2e-4", "3e-4"]
savenames = ["1e-8", "3e-8", "1e-7", "3e-7", "1e-6", "3e-6", "1e-5"]

leg = ["$d = 10$ nm", "$d = 30$ nm", "$d = 100$ nm", "$d = 300$ nm",
       "$d = 1$ \u03BCm", "$d = 3$ \u03BCm", "$d = 10$ \u03BCm"]
leg.append("rescaled")

sigm_0_3 = []
tau_0_3 = []
sigm_3 = []
tau_3 = []

sigm_0_2 = []
tau_0_2 = []
sigm_2 = []
tau_2 = []

specs_0_3 = []
specs_50_3 = []
specs_100_3 = []

for p in ["0", "50", "100"]:

    dat_3 = []
    dat_2 = []
    for filename in filenames:
        lspec_3 = (load_spec("Data/distance_dependence/parallel/p_"+p+"/3d",
                             "d_"+filename+".txt"))
        lspec_2 = (load_spec("Data/distance_dependence/parallel/p_"+p+"/2d",
                             "d_"+filename+".txt"))
        dat_3.append(lspec_3)
        dat_2.append(lspec_2)
        if p == "0" and filename in savenames:
            specs_0_3.append(lspec_3)
        elif p == "50" and filename in savenames:
            specs_50_3.append(lspec_3)
        elif p == "100" and filename in savenames:
            specs_100_3.append(lspec_3)
    
    dat_comp_3 = load_spec("Data/one_sphere/3d", "p_"+p+".txt")
    spec_rs_3 = scaleMaxwellGarnett(nu_3, dat_comp_3[1], nu0_3, dim="3d")
    
    dat_comp_2 = load_spec("Data/one_sphere/2d", "p_"+p+".txt")
    spec_rs_2 = scaleMaxwellGarnett(nu_2, dat_comp_2[1], nu0_2, dim="2d")
    
    if p == "0":
        specs_0_3.append((dat_comp_3[0], spec_rs_3))
    elif p == "50":
        specs_50_3.append((dat_comp_3[0], spec_rs_3))
    if p == "100":
        specs_100_3.append((dat_comp_3[0], spec_rs_3))
    
    tau_x_3 = []
    sigm_x_3 = []
    
    for spec in dat_3:
        w = spec[0]
        
        sig_temp, tau_temp = findmax(w, spec[1].imag)
        
        sigm_x_3.append(sig_temp)
        tau_x_3.append(tau_temp)
    
    sigm_os_3, tau_os_3 = findmax(w, spec_rs_3.imag)
    sigm_0_3.append(sigm_os_3)
    tau_0_3.append(tau_os_3)
        
    sigm_3.append(sigm_x_3)
    tau_3.append(tau_x_3)
    
    tau_x_2 = []
    sigm_x_2 = []
    
    for spec in dat_2:
        w = spec[0]
        
        sig_temp, tau_temp = findmax(w, spec[1].imag)
        
        sigm_x_2.append(sig_temp)
        tau_x_2.append(tau_temp)
    
    sigm_os_2, tau_os_2 = findmax(dat_comp_2[0], spec_rs_2.imag)
    sigm_0_2.append(sigm_os_2)
    tau_0_2.append(tau_os_2)
        
    sigm_2.append(sigm_x_2)
    tau_2.append(tau_x_2)
    

plt.rcParams.update({'font.size': 14})
plt.rcParams["figure.figsize"] = (14,9)

col = plt.cm.viridis(np.linspace(0.1, 1, 7))

fig, [[ax0, ax1], [ax2, ax3]] = plt.subplots(2,2)

ii = 0
for spec_0 in specs_0_3:
    if (spec_0[1] == specs_0_3[-1][1]).all():
        ax0.semilogx(spec_0[0], spec_0[1].real, "k--")
    else:
        ax0.semilogx(spec_0[0], spec_0[1].real, color=col[ii],
                     linewidth=2)
        ii += 1

ii = 0
for spec_0 in specs_0_3:
    if (spec_0[1] == specs_0_3[-1][1]).all():
        ax1.semilogx(spec_0[0], spec_0[1].imag, "k--")
    else:
        ax1.semilogx(spec_0[0], spec_0[1].imag, color=col[ii],
                     linewidth=2)
        ii += 1
ax1.set_ylim([-0.1e-7, 0.3e-6])

ii = 0
for spec_50 in specs_50_3:
    if (spec_50[1] == specs_50_3[-1][1]).all():
        ax2.semilogx(spec_50[0], spec_50[1].imag, "k--")
    else:
        ax2.semilogx(spec_50[0], spec_50[1].imag, color=col[ii],
                     linewidth=2)
        ii += 1

ii = 0
for spec_100 in specs_100_3:
    if (spec_100[1] == specs_100_3[-1][1]).all():
        ax3.semilogx(spec_100[0], spec_100[1].imag, "k--")
    else:
        ax3.semilogx(spec_100[0], spec_100[1].imag, color=col[ii],
                     linewidth=2)
        ii += 1
        
ax1.legend(leg)

ax0.set_xlabel("$\u03C9$ [rad/s]")
ax0.set_ylabel("$\u03C3'/\u03C3_0$ [-]")
ax1.set_xlabel("$\u03C9$ [rad/s]")
ax1.set_ylabel("$\u03C3''/\u03C3_0$ [-]")
ax2.set_xlabel("$\u03C9$ [rad/s]")
ax2.set_ylabel("$\u03C3$''/$\u03C3_0$ [-]")
ax3.set_xlabel("$\u03C9$ [rad/s]")
ax3.set_ylabel("$\u03C3''/\u03C3_0$ [-]")

ax0.set_title("(a)", loc="left")
ax1.set_title("(b)", loc="left")
ax2.set_title("(c)", loc="left")
ax3.set_title("(d)", loc="left")

fig.tight_layout()
fig.savefig("Figures/para_spec.png", dpi=300, bbox_inches="tight")

plt.rcParams["figure.figsize"] = (14,9)

d = np.array([float(name) for name in filenames])

fig, ax = plt.subplots(2,2)
ax[0,1].semilogx(d, sigm_3[0], "o-", color=col[6])
ax[0,1].semilogx(d, sigm_3[1], "o-", color=col[3])
ax[0,1].semilogx(d, sigm_3[2], "o-", color=col[0])
ax[0,1].semilogx([d[0], d[-1]], [sigm_0_3[0], sigm_0_3[0]], "--",
                 color=col[6])
ax[0,1].semilogx([d[0], d[-1]], [sigm_0_3[1], sigm_0_3[1]], "--",
                 color=col[3])
ax[0,1].semilogx([d[0], d[-1]], [sigm_0_3[2], sigm_0_3[2]], "--",
                 color=col[0])

ax[1,1].semilogx(d, sigm_2[0], "o-", color=col[6])
ax[1,1].semilogx(d, sigm_2[1], "o-", color=col[3])
ax[1,1].semilogx(d, sigm_2[2], "o-", color=col[0])
ax[1,1].semilogx([d[0], d[-1]], [sigm_0_2[0], sigm_0_2[0]], "--",
                 color=col[6])
ax[1,1].semilogx([d[0], d[-1]], [sigm_0_2[1], sigm_0_2[1]], "--",
                 color=col[3])
ax[1,1].semilogx([d[0], d[-1]], [sigm_0_2[2], sigm_0_2[2]], "--",
                 color=col[0])

ax[0,0].loglog(d, tau_3[0], "-o", color=col[6])
ax[0,0].loglog(d, tau_3[1], "-o", color=col[3])
ax[0,0].loglog(d, tau_3[2], "-o", color=col[0])
ax[0,0].loglog([d[0], d[-1]], [tau_0_3[0], tau_0_3[0]], "--", color=col[6])
ax[0,0].loglog([d[0], d[-1]], [tau_0_3[1], tau_0_3[1]], "--", color=col[3])
ax[0,0].loglog([d[0], d[-1]], [tau_0_3[2], tau_0_3[2]], "--", color=col[0])

ax[1,0].loglog(d, tau_2[0], "-o", color=col[6])
ax[1,0].loglog(d, tau_2[1], "-o", color=col[3])
ax[1,0].loglog(d, tau_2[2], "-o", color=col[0])
ax[1,0].loglog([d[0], d[-1]], [tau_0_2[0], tau_0_2[0]], "--", color=col[6])
ax[1,0].loglog([d[0], d[-1]], [tau_0_2[1], tau_0_2[1]], "--", color=col[3])
ax[1,0].loglog([d[0], d[-1]], [tau_0_2[2], tau_0_2[2]], "--", color=col[0])

ax[1,0].legend(["$p = 0.0$", "$p = 0.5$", "$p = 1.0$"], loc='lower right')

ax[0,1].set_xlabel("$d$ [m]")
ax[0,1].set_ylabel("$\u03C3_{max}''/\u03C3_0$ [-]")
ax[1,1].set_xlabel("$d$ [m]")
ax[1,1].set_ylabel("$\u03C3_{max}''/\u03C3_0$ [-]")
ax[0,0].set_xlabel("$d$ [m]")
ax[0,0].set_ylabel("$\u03C4$ [s]")
ax[1,0].set_xlabel("$d$ [m]")
ax[1,0].set_ylabel("$\u03C4$ [s]")

ax[0,0].set_title("(a)", loc="left")
ax[0,1].set_title("(b)", loc="left")
ax[1,0].set_title("(c)", loc="left")
ax[1,1].set_title("(d)", loc="left")

fig.tight_layout()
fig.savefig("Figures/para_comp.png",  dpi=300, bbox_inches="tight")
