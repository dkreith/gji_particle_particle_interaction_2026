import numpy as np
import matplotlib.pyplot as plt
from functions import load_spec

fnames = ["3", "4", "5", "6"]

datwth = []
datwot = []
for fname in fnames:
    datwth.append(load_spec("Data/finite_cluster/overlap",
                            "size_"+fname+"_with.txt"))
    datwot.append(load_spec("Data/finite_cluster/overlap",
                            "size_"+fname+"_without.txt"))
    
col = plt.cm.viridis(np.linspace(0.1, 1, 4))

fig, ax = plt.subplots(figsize=(8,6.5))
ax.loglog(2*np.pi*datwth[0][0], datwth[0][1].imag, color=col[0], linewidth=2)
ax.loglog(2*np.pi*datwot[0][0], datwot[0][1].imag, "--", color=col[0],
          linewidth=2)
ax.loglog(2*np.pi*datwth[1][0], datwth[1][1].imag, color=col[1], linewidth=2)
ax.loglog(2*np.pi*datwot[1][0], datwot[1][1].imag, "--", color=col[1],
          linewidth=2)
ax.loglog(2*np.pi*datwth[2][0], datwth[2][1].imag, color=col[2], linewidth=2)
ax.loglog(2*np.pi*datwot[2][0], datwot[2][1].imag, "--", color=col[2],
          linewidth=2)
ax.loglog(2*np.pi*datwth[3][0], datwth[3][1].imag, color=col[3], linewidth=2)
ax.loglog(2*np.pi*datwot[3][0], datwot[3][1].imag, "--", color=col[3],
          linewidth=2)
ax.legend(["3x3 cluster", "3x3 eq. particle",
           "4x4 cluster", "4x4 eq. particle",
           "5x5 cluster", "5x5 eq. particle",
           "6x6 cluster", "6x6 eq. particle",])
ax.set_xlabel("$\u03C9$ [rad/s]")
ax.set_ylabel("$\u03C3''/\u03C3_0$ [-]")
fig.tight_layout()
fig.savefig("Figures/finite_cluster_overlap.pdf", dpi=300, bbox_inches="tight")