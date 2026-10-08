import numpy as np
import matplotlib.pyplot as plt
from functions import load_spec

fnames = ["1,99", "1,75", "1,50", "1,25", "1,00", "0,75", "0,50", "0,25"]

dat_para = []
dat_perp = []
for fname in fnames:
    dat_para.append(load_spec("Data/overlap_two_spheres/parallel",
                              "dc_"+fname+"_a.txt"))
    dat_perp.append(load_spec("Data/overlap_two_spheres/perpendicular",
                              "dc_"+fname+"_a.txt"))

dat0 = load_spec("Data/one_sphere/2d/", "p_50.txt")
    
col = plt.cm.viridis(np.linspace(0.1, 0.95, len(fnames)))
plt.rcParams.update({'font.size': 14})
plt.rcParams["figure.figsize"] = (9,8)
    
fig, ax = plt.subplots(2,1)
for ii in range(len(fnames)):
    ax[0].semilogx(dat_para[ii][0], dat_para[ii][1].imag, color=col[ii],
                   linewidth=2)
    ax[1].semilogx(dat_perp[ii][0], dat_perp[ii][1].imag, color=col[ii],
                   linewidth=2)
ax[0].semilogx(dat0[0], dat0[1].imag, "--k", linewidth=2)
ax[1].semilogx(dat0[0], dat0[1].imag, "--k", linewidth=2)
ax[0].set_ylim([0, 0.000017])
ax[1].set_ylim([0, 0.000017])
ax[1].legend(["$d_c = 1.99 a$", "$d_c = 1.75 a$", "$d_c = 1.5 a$",
              "$d_c = 1.25 a$", "$d_c = a$", "$d_c = 0.75 a$", "$d_c = 0.5 a$",
              "$d_c = 0.25 a$", "$d_c = 0$"], loc="upper right", ncol=2)
ax[0].set_ylabel("$\u03C3''/\u03C3_0$ [-]")
ax[1].set_ylabel("$\u03C3''/\u03C3_0$ [-]")
ax[0].set_xlabel("$\u03C9$ [rad/s]")
ax[1].set_xlabel("$\u03C9$ [rad/s]")
ax[0].set_title("(a)", loc="left")
ax[1].set_title("(b)", loc="left")
fig.tight_layout()
fig.savefig("Figures/overlap_two_spheres.pdf", dpi=300)