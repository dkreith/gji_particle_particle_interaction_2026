import numpy as np
import matplotlib.pyplot as plt
from functions import load_spec, findmax, scaleMaxwellGarnett

Vs = np.pi*(5e-6)**2
Vd = (5e-5)**2
nu = Vs/Vd
nu0 = 2*nu

filenames = ["1e-8", "15e-9", "2e-8", "3e-8", "5e-8", "7e-8",
             "1e-7", "15e-8", "2e-7", "3e-7", "5e-7", "7e-7",
             "1e-6", "15e-7", "2e-6", "3e-6", "5e-6", "7e-6",
             "1e-5", "15e-6", "2e-5"]
savenames = ["1e-8", "3e-8", "1e-7", "3e-7", "1e-6", "3e-6", "1e-5"]

leg = ["$d = 10$ nm", "$d = 30$ nm", "$d = 100$ nm", "$d = 300$ nm",
       "$d = 1$ \u03BCm", "$d = 3$ \u03BCm", "$d = 10$ \u03BCm"]
leg.append("rescaled")

sigm_0 = []
sigm = []
tau_0 = []
tau = []

specs_0 = []
specs_50 = []
specs_100 = []

for p in ["0", "50", "100"]:

    dat = []
    for filename in filenames:
        lspec = (load_spec("Data/distance_dependence/perpendicular/p_"+p,
                            "d_"+filename+".txt"))
        dat.append(lspec)
        if p == "0" and filename in savenames:
            specs_0.append(lspec)
        elif p == "50" and filename in savenames:
            specs_50.append(lspec)
        elif p == "100" and filename in savenames:
            specs_100.append(lspec)
    
    dat_comp = load_spec("Data/one_sphere/2d", "p_"+p+".txt")
    
    spec_rs = scaleMaxwellGarnett(nu, dat_comp[1], nu0, dim="3d")
    
    if p == "0":
        specs_0.append((dat_comp[0], spec_rs))
    elif p == "50":
        specs_50.append((dat_comp[0], spec_rs))
    if p == "100":
        specs_100.append((dat_comp[0], spec_rs))
    
    tau_x = []
    sigm_x = []
    
    for spec in dat:
        w = spec[0]
        
        sig_temp, tau_temp = findmax(w, spec[1].imag)
        
        sigm_x.append(sig_temp)
        tau_x.append(tau_temp)
    
    sigm_os, tau_os = findmax(w, spec_rs.imag)
    sigm_0.append(sigm_os)
    tau_0.append(tau_os)
        
    sigm.append(sigm_x)
    tau.append(tau_x)
    
plt.rcParams.update({'font.size': 14})
plt.rcParams["figure.figsize"] = (13,13)

col = plt.cm.viridis(np.linspace(0.1, 1, 7))

d = np.array([float(name) for name in filenames])

fig, ax = plt.subplot_mosaic([["spec0", "sigm"],
                              ["spec0", "sigm"],
                              ["spec50", "sigm"],
                              ["spec50", "tau"],
                              ["spec100", "tau"],
                              ["spec100", "tau"]])

ii = 0
for spec_0 in specs_0:
    if (spec_0[1] == specs_0[-1][1]).all():
        ax["spec0"].semilogx(spec_0[0], spec_0[1].imag*1e5, "k--")
    else:
        ax["spec0"].semilogx(spec_0[0], spec_0[1].imag*1e5, color=col[ii],
                         linewidth=2)
        ii += 1

ii = 0
for spec_50 in specs_50:
    if (spec_50[1] == specs_50[-1][1]).all():
        ax["spec50"].semilogx(spec_50[0], spec_50[1].imag*1e5, "k--")
    else:
        ax["spec50"].semilogx(spec_50[0], spec_50[1].imag*1e5, color=col[ii],
                         linewidth=2)
        ii += 1

ii = 0
for spec_100 in specs_100:
    if (spec_100[1] == specs_100[-1][1]).all():
        ax["spec100"].semilogx(spec_100[0], spec_100[1].imag*1e5, "k--")
    else:
        ax["spec100"].semilogx(spec_100[0], spec_100[1].imag*1e5, color=col[ii],
                         linewidth=2)
        ii += 1
        
ax["spec0"].legend(leg)

ax["sigm"].semilogx(d, np.array(sigm[0])*1e5, "-o", color=col[6])
ax["sigm"].semilogx(d[5:-1], np.array(sigm[1][5:-1])*1e5, "-o", color=col[3])
ax["sigm"].semilogx(d[3:-1], np.array(sigm[2][3:-1])*1e5, "-o", color=col[0])
ax["sigm"].semilogx(d[0:6], np.array(sigm[1][0:6])*1e5,
                    linestyle=(0, (3,1,1,1)), color=col[3])
ax["sigm"].semilogx(d[0:6], np.array(sigm[1][0:6])*1e5, "o", color=col[3])
ax["sigm"].semilogx(d[0:4], np.array(sigm[2][0:4])*1e5,
                    linestyle=(0, (3,1,1,1)), color=col[0])
ax["sigm"].semilogx(d[0:4], np.array(sigm[2][0:4])*1e5, "o", color=col[0])
ax["sigm"].semilogx([d[0], d[-1]], [sigm[0][-1]*1e5, sigm[0][-1]*1e5], "--",
                 color=col[6])
ax["sigm"].semilogx([d[0], d[-1]], [sigm[1][-1]*1e5, sigm[1][-1]*1e5], "--",
                 color=col[3])
ax["sigm"].semilogx([d[0], d[-1]], [sigm[2][-1]*1e5, sigm[2][-1]*1e5], "--",
                 color=col[0])
ax["sigm"].legend(["$p = 0.0$", "$p = 0.5$", "$p = 1.0$"])

ax["tau"].loglog(d, tau[0], "-o", color=col[6])
ax["tau"].loglog(d[5:-1], tau[1][5:-1], "-o", color=col[3])
ax["tau"].loglog(d[3:-1], tau[2][3:-1], "-o", color=col[0])
ax["tau"].loglog(d[0:5], tau[1][0:5], linestyle=(0, (3,1,1,1)), color=col[3])
ax["tau"].loglog(d[0:5], tau[1][0:5], "o", color=col[3])
ax["tau"].loglog(d[0:3], tau[2][0:3], linestyle=(0, (3,1,1,1)), color=col[0])
ax["tau"].loglog(d[0:3], tau[2][0:3], "o", color=col[0])
ax["tau"].loglog([d[0], d[-1]], [tau_0[0], tau_0[0]], "--", color=col[6])
ax["tau"].loglog([d[0], d[-1]], [tau_0[1], tau_0[1]], "--", color=col[3])
ax["tau"].loglog([d[0], d[-1]], [tau_0[2], tau_0[2]], "--", color=col[0])

ax["spec0"].set_xlabel("$\u03C9$ [rad/s]")
ax["spec0"].set_ylabel("$\u03C3''/\u03C3_0$ [-]")
ax["spec50"].set_xlabel("$\u03C9$ [rad/s]")
ax["spec50"].set_ylabel("$\u03C3$''/$\u03C3_0$ [-]")
ax["spec100"].set_xlabel("$\u03C9$ [rad/s]")
ax["spec100"].set_ylabel("$\u03C3''/\u03C3_0$ [-]")
ax["sigm"].set_xlabel("$d$ [m]")
ax["sigm"].set_ylabel("$\u03C3_{max}''/\u03C3_0$ [-]")
ax["tau"].set_xlabel("$d$ [m]")
ax["tau"].set_ylabel("$\u03C4$ [s]")

ax["spec0"].set_title("(a)\n1e-5", loc="left")
ax["spec50"].set_title("(b)\n1e-5", loc="left")
ax["spec100"].set_title("(c)\n1e-5", loc="left")
ax["sigm"].set_title("(d)\n1e-5", loc="left")
ax["tau"].set_title("(e)", loc="left")

fig.tight_layout()
fig.savefig("Figures/perp.png",  dpi=300, bbox_inches="tight")