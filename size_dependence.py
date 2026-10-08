import numpy as np
import matplotlib.pyplot as plt
from functions import load_spec

a = ["5e-7", "1e-6", "2,5e-6", "5e-6", "1e-5", "2,5e-5", "5e-5"]
d = ["1e-8", "2e-8", "5e-8", "1e-7", "2e-7", "5e-7", "1e-6"]

a_a = [5e-7, 1e-6, 2.5e-6, 5e-6, 1e-5, 2.5e-5, 5e-5]

dat_rat = []
dat_gap = []

for ii in range(len(a)):
    dat_rat.append(load_spec("Data/size_dependence/same_ratio",
                             "a_"+a[ii]+"_d_"+d[ii]+"_phi_90.txt"))
    dat_gap.append(load_spec("Data/size_dependence/same_gap",
                             "a_"+a[ii]+"_d_1e-7_phi_90.txt"))
    
leg = ["$a = 0.5$ \u03BCm", "$a = 1$ \u03BCm", "$a = 2.5$ \u03BCm",
       "$a = 5$ \u03BCm", "$a = 10$ \u03BCm", "$a = 25$ \u03BCm",
       "$a = 50$ \u03BCm"]

plt.rcParams["figure.figsize"] = (9,9.5)

col = plt.cm.viridis(np.linspace(0.1, 1, 7))

fig, ax = plt.subplots(2,1)
for ii in range(len(a_a)):
    ax[0].semilogx(dat_rat[ii][0], dat_rat[ii][1].imag*a_a[ii]/5e-6*1e5,
                   linewidth=2, color=col[ii])
    ax[1].semilogx(dat_gap[ii][0], dat_gap[ii][1].imag*a_a[ii]/5e-6*1e5,
                   linewidth=2, color=col[ii])
ax[0].set_xlabel("$\u03C9$ [rad/s]")
ax[0].set_ylabel("$\u03C3''/\u03C3_0 \cdot a/a_0$ [-]")
ax[1].set_xlabel("$\u03C9$ [rad/s]")
ax[1].set_ylabel("$\u03C3''/\u03C3_0 \cdot a/a_0$ [-]")
ax[0].set_ylim([-5e-1, 10])
ax[1].set_ylim([-1e-1, 20])
ax[0].set_title("(a)\n1e-5", loc="left")
ax[1].set_title("(b)\n1e-5", loc="left")
ax[1].legend(leg, loc="upper right")
fig.tight_layout()

fig.savefig("Figures/size_dependence.pdf", dpi=300, bbox_inches="tight")