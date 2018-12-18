import numpy as np
import matplotlib.pyplot as plt
plotphi1 = np.array([2.68915606, 2.68915606,,3.1795165,  1.52639784, 0.35340126, 0.20909307 ,0.48643007,
  0.21986978, 0.22957882 ,0.51787338 ,0.62186825])
plotphi2 = np.array([2.44736662,2.44736662,3.04602345, 1.55553508 ,0.37605368 ,0.22908949, 0.49570432,
  0.23040307 ,0.24666067 ,0.57582504 ,0.76415354])
plotphi1 = np.flipud(plotphi1)
plotphi2 = np.flipud(plotphi2)
# group_grid = np.log([0, 0.058,0.14,0.28,0.625,4,10,40,5530,821e3,20e6])
group_grid = [0.0001, 0.058,0.14,0.28,0.625,4,10,40,5530,821e3,20e6]
plt.step( group_grid[0:11], plotphi1,where='post', label='fuel')
plt.step(group_grid[0:11], plotphi2, where='post', label='mod')
plt.legend()
plt.xlim(0.0001,20e6)
# plt.ylim(5,15)
plt.xlabel('Energy')
plt.ylabel('Scalar Flux')
plt.xscale('log')
plt.show()