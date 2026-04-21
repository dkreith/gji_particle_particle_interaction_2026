import numpy as np
import matplotlib.pyplot as plt
from functions import load_spec, findmax, scaleMaxwellGarnett

def sine_fun(phi0, sigm):
    phi = phi0 * np.pi/180
    A = sigm[0] - sigm[-1]
    b = sigm[-1]
    f = A*(np.cos(2*phi) + 1)/2 + b
    return f

"""
distance: 1e-6
sphere radius: 5e-6
sig: 0.01
"""

Vs = np.pi*(5e-6)**2
Vd = 5e-5**2
nu = Vs/Vd
nu0 = 2*nu

leg = ["$\u03B1 = 0$°", "$\u03B1 = 15$°", "$\u03B1 = 30$°", "$\u03B1 = 45$°",
       "$\u03B1 = 60$°", "$\u03B1 = 75$°", "$\u03B1 = 90$°"]
leg.append("rescaled")

filenames = ["0", "5", "10", "15", "20", "25", "30", "35", "40", "45",
             "50", "55", "60", "65", "70", "75", "80", "85", "90"]
savenames = ["0", "15", "30", "45", "60", "75", "90"]

dnames = ["1e-7", "3e-7", "1e-6", "3e-6"]

sigm_00 = []
sigm_50 = []
sigm_100 = []

plt.rcParams.update({'font.size': 14})
plt.rcParams["figure.figsize"] = (13,9)

fig0, ax0 = plt.subplots(2,2)

for dn in dnames:
    
    sigm_0 = []
    sigm = []
    tau = []
    
    specs_0 = []
    specs_50 = []
    specs_100 = []
    
    for p in ["0", "50", "100"]:
    
        dat = []
        for filename in filenames:
            lspec = (load_spec("Data/angle_dependence/d_"+dn+"/two_spheres_p_"+p,
                                "d_"+dn+"_phi_"+filename+".txt"))
            dat.append(lspec)
            if p == "0" and filename in savenames:
                specs_0.append(lspec)
            elif p == "50" and filename in savenames:
                specs_50.append(lspec)
            elif p == "100" and filename in savenames:
                specs_100.append(lspec)
            
        
        dat_comp = load_spec("Data/one_sphere/2d", "p_"+p+".txt")
        
        spec_rs = scaleMaxwellGarnett(nu, dat_comp[1], nu0, dim="2d")
        
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
        
        sigm.append(sigm_x)
        tau.append(tau_x)
    
    sigm_00.append(sigm[0])
    sigm_50.append(sigm[1])
    sigm_100.append(sigm[2])
    
    phi = np.array([float(name) for name in filenames])
    
    col = plt.cm.viridis(np.linspace(0.1, 1, 7))
    
    fig, ax = plt.subplots(2,2)
    
    ii = 6
    for spec_0 in specs_0:
        if (spec_0[1] == specs_0[-1][1]).all():
            ax[0,0].semilogx(spec_0[0], spec_0[1].imag, "k--")
        else:
            ax[0,0].semilogx(spec_0[0], spec_0[1].imag, color=col[ii],
                             linewidth=2)
            ii -= 1
    
    if dn == "1e-6":
        ax[0,0].set_ylim([-0.1e-6, 1.1e-5])
    elif dn == "1e-7":
        ax[0,0].set_ylim([-0.1e-6, 1e-4])
    elif dn == "3e-7":
        ax[0,0].set_ylim([-0.1e-6, 4e-5])
    elif dn == "3e-6":
        ax[0,0].set_ylim([-0.1e-6, 5e-6])
        
    ii = 6
    for spec_50 in specs_50:
        if (spec_50[1] == specs_50[-1][1]).all():
            ax[0,1].semilogx(spec_50[0], spec_50[1].imag, "k--")
        else:
            ax[0,1].semilogx(spec_50[0], spec_50[1].imag, color=col[ii],
                             linewidth=2)
            ii -= 1
    
    if dn == "1e-6":
        ii = 6
        for spec_0 in specs_0:
            if (spec_0[1] == specs_0[-1][1]).all():
                ax0[0,0].semilogx(spec_0[0], spec_0[1].imag, "k--")
            else:
                ax0[0,0].semilogx(spec_0[0], spec_0[1].imag, color=col[ii],
                                  linewidth=2)
                ii -= 1
    elif dn == "1e-7":
        ii = 6
        for spec_0 in specs_0:
            if (spec_0[1] == specs_0[-1][1]).all():
                ax0[1,0].semilogx(spec_0[0], spec_0[1].imag, "k--")
            else:
                ax0[1,0].semilogx(spec_0[0], spec_0[1].imag, color=col[ii],
                                  linewidth=2)
                ii -= 1
    if dn == "1e-6":
        ii = 6
        for spec_50 in specs_50:
            if (spec_50[1] == specs_50[-1][1]).all():
                ax0[0,1].semilogx(spec_50[0], spec_50[1].imag, "k--")
            else:
                ax0[0,1].semilogx(spec_50[0], spec_50[1].imag, color=col[ii],
                                  linewidth=2)
                ii -= 1
    elif dn == "1e-7":
        ii = 6
        for spec_50 in specs_50:
            if (spec_50[1] == specs_50[-1][1]).all():
                ax0[1,1].semilogx(spec_50[0], spec_50[1].imag, "k--")
            else:
                ax0[1,1].semilogx(spec_50[0], spec_50[1].imag, color=col[ii],
                                  linewidth=2)
                ii -= 1
    
    ii = 6
    for spec_100 in specs_100:
        if (spec_100[1] == specs_100[-1][1]).all():
            ax[1,0].semilogx(spec_100[0], spec_100[1].imag, "k--")
        else:
            ax[1,0].semilogx(spec_100[0], spec_100[1].imag, color=col[ii],
                             linewidth=2)
            ii -= 1
    
    if dn == "1e-6":
        ax[0,1].legend(leg, loc="upper right")
        ax0[0,0].legend(leg, loc="upper right")
    elif dn == "1e-7":
        ax[0,0].legend(leg, loc="upper left")
    elif dn == "3e-7":
        ax[0,0].legend(leg, loc="upper left")
    elif dn == "3e-6":
        ax[0,1].legend(leg, loc="upper right")
    
    ax[1,1].plot(phi, sigm[0], "o-", color=col[6], linewidth=2)
    ax[1,1].plot(phi, sigm[1], "o-", color=col[3], linewidth=2)
    ax[1,1].plot(phi, sigm[2], "o-", color=col[0], linewidth=2)
    ax[1,1].plot(phi, sine_fun(phi, sigm[0]), color=col[6], linewidth=2)
    ax[1,1].plot(phi, sine_fun(phi, sigm[1]), color=col[3], linewidth=2)
    ax[1,1].plot(phi, sine_fun(phi, sigm[2]), color=col[0], linewidth=2)
    ax[1,1].legend(["$p = 0.0$", "$p = 0.5$", "$p = 1.0$"])
    ax[1,1].plot([phi[0], phi[-1]], [sigm_0[0], sigm_0[0]], "--",
                 color=col[6])
    ax[1,1].plot([phi[0], phi[-1]], [sigm_0[1], sigm_0[1]], "--",
                 color=col[3])
    ax[1,1].plot([phi[0], phi[-1]], [sigm_0[2], sigm_0[2]], "--",
                 color=col[0])
    
    ax[0,0].set_xlabel("$\u03C9$ [rad/s]")
    ax[0,0].set_ylabel("$\u03C3''/\u03C3_0$ [-]")
    ax[1,0].set_xlabel("$\u03C9$ [rad/s]")
    ax[1,0].set_ylabel("$\u03C3$''/$\u03C3_0$ [-]")
    ax[0,1].set_xlabel("$\u03C9$ [rad/s]")
    ax[0,1].set_ylabel("$\u03C3''/\u03C3_0$ [-]")
    ax[1,1].set_xlabel("$\u03B1$ [°]")
    ax[1,1].set_ylabel("$\u03C3_{max}''/\u03C3_0$ [-]")
    
    ax[0,0].set_title("(a)", loc="left")
    ax[0,1].set_title("(b)", loc="left")
    ax[1,0].set_title("(c)", loc="left")
    ax[1,1].set_title("(d)", loc="left")
    
    fig.tight_layout()
    fig.savefig("Figures/angle_"+dn+".png", dpi=300)

ax0[0,0].set_ylim([-0.1e-6, 1.1e-5])

ax0[0,0].set_xlabel("$\u03C9$ [rad/s]")
ax0[0,0].set_ylabel("$\u03C3''/\u03C3_0$ [-]")
ax0[1,0].set_xlabel("$\u03C9$ [rad/s]")
ax0[1,0].set_ylabel("$\u03C3''/\u03C3_0$ [-]")
ax0[0,1].set_xlabel("$\u03C9$ [rad/s]")
ax0[0,1].set_ylabel("$\u03C3''/\u03C3_0$ [-]")
ax0[1,1].set_xlabel("$\u03C9$ [rad/s]")
ax0[1,1].set_ylabel("$\u03C3''/\u03C3_0$ [-]")

ax0[0,0].set_title("(a)", loc="left")
ax0[0,1].set_title("(b)", loc="left")
ax0[1,0].set_title("(c)", loc="left")
ax0[1,1].set_title("(d)", loc="left")    

fig0.tight_layout()
fig0.savefig("Figures/angle_spec.png", dpi=300)
   
plt.rcParams["figure.figsize"] = (6,11)

fig, ax = plt.subplots(3,1)
ax[0].plot(phi, sigm_00[0], "o-", color=col[0], linewidth=2)
ax[0].plot(phi, sigm_00[1], "o-", color=col[2], linewidth=2)
ax[0].plot(phi, sigm_00[2], "o-", color=col[4], linewidth=2)
ax[0].plot(phi, sigm_00[3], "o-", color=col[6], linewidth=2)
ax[0].plot([phi[0], phi[-1]], [sigm_0[0], sigm_0[0]], "--", color="gray")
ax[0].plot(phi, sine_fun(phi, sigm_00[0]), color=col[0], linewidth=2)
ax[0].plot(phi, sine_fun(phi, sigm_00[1]), color=col[2], linewidth=2)
ax[0].plot(phi, sine_fun(phi, sigm_00[2]), color=col[4], linewidth=2)
ax[0].plot(phi, sine_fun(phi, sigm_00[3]), color=col[6], linewidth=2)
ax[1].plot(phi, sigm_50[0], "o", color=col[0], linewidth=2)
ax[1].plot(phi, sigm_50[1], "o", color=col[2], linewidth=2)
ax[1].plot(phi, sigm_50[2], "o", color=col[4], linewidth=2)
ax[1].plot(phi, sigm_50[3], "o", color=col[6], linewidth=2)
ax[1].plot([phi[0], phi[-1]], [sigm_0[1], sigm_0[1]], "--", color="gray")
ax[1].plot(phi, sine_fun(phi, sigm_50[0]), color=col[0], linewidth=2)
ax[1].plot(phi, sine_fun(phi, sigm_50[1]), color=col[2], linewidth=2)
ax[1].plot(phi, sine_fun(phi, sigm_50[2]), color=col[4], linewidth=2)
ax[1].plot(phi, sine_fun(phi, sigm_50[3]), color=col[6], linewidth=2)
ax[2].plot(phi, sigm_100[0], "o", color=col[0], linewidth=2)
ax[2].plot(phi, sigm_100[1], "o", color=col[2], linewidth=2)
ax[2].plot(phi, sigm_100[2], "o", color=col[4], linewidth=2)
ax[2].plot(phi, sigm_100[3], "o", color=col[6], linewidth=2)
ax[2].plot([phi[0], phi[-1]], [sigm_0[2], sigm_0[2]], "--", color="gray")
ax[2].plot(phi, sine_fun(phi, sigm_100[0]), color=col[0], linewidth=2)
ax[2].plot(phi, sine_fun(phi, sigm_100[1]), color=col[2], linewidth=2)
ax[2].plot(phi, sine_fun(phi, sigm_100[2]), color=col[4], linewidth=2)
ax[2].plot(phi, sine_fun(phi, sigm_100[3]), color=col[6], linewidth=2)

ax[0].set_xlabel("$\u03B1$ [°]")
ax[0].set_ylabel("$\u03C3_{max}''/\u03C3_0$ [-]")
ax[1].set_xlabel("$\u03B1$ [°]")
ax[1].set_ylabel("$\u03C3_{max}''/\u03C3_0$ [-]")
ax[2].set_xlabel("$\u03B1$ [°]")
ax[2].set_ylabel("$\u03C3_{max}''/\u03C3_0$ [-]")

ax[0].set_title("(a)", loc="left")
ax[1].set_title("(b)", loc="left")
ax[2].set_title("(c)", loc="left")

ax[0].legend(["$d = 0.1$ \u03BCm", "$d = 0.3$ \u03BCm", "$d = 1$ \u03BCm",
              "$d = 3$ \u03BCm", "rescaled"], loc="upper left")

fig.tight_layout()
fig.savefig("Figures/angle_max.png", dpi=300)