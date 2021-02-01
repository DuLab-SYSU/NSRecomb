import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mtp

plt.style.use('ggplot')

mu = np.array([0.25, 1.25, 2.5, 5])
mu_1 = np.array([1, 2, 3, 4])
r = np.array([0, 0.25, 1, 4, 16])
r_1 = np.array([0, 1 ,2 ,3 ,4])

Z = np.array([[0, 0, 2, 5],
              [0, 4, 9, 14],
              [0, 12, 23, 29],
              [0, 32, 58, 77],
              [1, 59, 91, 98]])


X, Y = np.meshgrid(mu, r)


fig = plt.figure(figsize=(7, 3.5))
ax = fig.add_subplot(121)

CS = ax.contourf(X, Y, Z, 20)
cbar = fig.colorbar(CS)
# ax.set_xlim(0.25, 5)
ax.set_xticks([0.25, 1.25, 2.5, 5])
# ax.set_ylim(0, 16)
ax.set_yticks([0.25, 1, 4, 16])
ax.set_xlabel("Mutation Rate ($10^{-5}$)", fontsize=16, font='Arial', color='black')
ax.set_ylabel("Reconbiantion Rate ($10^{-6}$)", fontsize=16, font='Arial', color='black')
# ax.set_title(interp)


ax = fig.add_subplot(122)

# for i, z_i in enumerate(Z):
#     ax.plot(mu, z_i, label="$r$=" + str(r[i]), linewidth=3)
#     ax.scatter(x=mu, y=z_i, s=50)

ax.plot(mu, [0, 0, 2, 5], linewidth=3, linestyle=':', color='C0', label='No recombination, our method')
ax.scatter(x=mu, y=[0, 0, 2, 5], color='C0', marker='s', s=50)

ax.plot(mu, [0, 2, 3, 23], linewidth=3, linestyle='-', color='C0', label='Convergent, our method')
ax.scatter(x=mu, y=[0, 2, 3, 23], color='C0', marker='s', s=50)

ax.plot(mu, [0, 0, 0, 1], linewidth=3, linestyle=':', color='C1', label='No recombination, 3seq')
ax.scatter(x=mu, y=[0, 0, 0, 1], color='C1', marker='s', s=50)

ax.plot(mu, [0, 3, 19, 35], linewidth=3, linestyle='-', color='C1', label='Convergent, 3seq')
ax.scatter(x=mu, y=[0, 3, 19, 35], color='C1', marker='s', s=50)

ax.set_xticks([0.25, 1.25, 2.5, 5])
ax.set_yticks([0, 20, 40, 60, 80, 100])
ax.set_xlabel("Mutation Rate ($10^{-5}$)", fontsize=16, font='Arial', color='black')
ax.set_ylabel("Probability of detecting\nrecombination", fontsize=16, font='Arial', color='black')
ax.legend()
# ax.legend(loc=2, bbox_to_anchor=(1, 1))

plt.show()