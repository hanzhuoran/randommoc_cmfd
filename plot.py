import numpy as np
import matplotlib.pyplot as plt
from util import *

# k = 1.10568
phidf = np.array([0.0499228759934832,0.0941120923671741,0.136239071905262,0.175380726724590,\
    0.210679405146194,0.241361683844689,0.266755309974487,0.286303922970518,0.299579235057431,0.306290405556839])
# k = 1.06048 
phirrF = np.array([0.27245382,0.48566274,0.68252428,0.86402302,1.0236056,1.16100592,\
    1.27493175,1.36304535,1.42169127,1.4510881])
# k = 1.08461
phirrL = np.array([0.27221076,0.48522969,0.68191615,0.86325392,1.02269544,1.15997473,\
 1.27380052,1.36183697,1.42043166,1.44980284])

phidf = phidf/np.linalg.norm(phidf)
phirrF = phirrF/np.linalg.norm(phirrF)
phirrL = phirrL/np.linalg.norm(phirrL)
axis = 1/2*(xbd_list[0:-1]+xbd_list[1:NX+1])
print(axis)
plt.plot(axis,phidf,label="diffusion")
plt.plot(axis,phirrF,'rx-',label="Fission iteration")
plt.plot(axis,phirrL,'bo-',label="Leakage iteration")
plt.legend()
plt.show()
# plotphi1 = np.flipud(plotphi1)
# plotphi2 = np.flipud(plotphi2)
# group_grid = np.log([0, 0.058,0.14,0.28,0.625,4,10,40,5530,821e3,20e6])
# group_grid = [0.0001, 0.058,0.14,0.28,0.625,4,10,40,5530,821e3,20e6]
# plt.step( group_grid[0:11], plotphi1,where='post', label='fuel')
# plt.step(group_grid[0:11], plotphi2, where='post', label='mod')
# plt.legend()
# plt.xlim(0.0001,20e6)
# plt.ylim(5,15)
# plt.xlabel('Energy')
# plt.ylabel('Scalar Flux')
# plt.xscale('log')
# plt.show()