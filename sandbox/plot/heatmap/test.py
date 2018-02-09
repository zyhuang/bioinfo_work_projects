import matplotlib.pyplot as plt
import numpy as np
from string import ascii_lowercase

fig, ax = plt.subplots()

a = np.random.random((16, 16))
cplot = ax.imshow(a, cmap='hot_r', interpolation='nearest')
cbar = fig.colorbar(cplot)

label = list(ascii_lowercase)[:16]
ax.set_xticks(np.arange(len(label)))
ax.set_yticks(np.arange(len(label)))
ax.set_xticklabels(label, minor=False)
ax.set_yticklabels(label, minor=False)

plt.show()
