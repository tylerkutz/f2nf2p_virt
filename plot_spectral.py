import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
from scipy.interpolate import griddata
from matplotlib import ticker, cm



sp_p = [];
alpha = [];
nu = [];
i = 0 
with open("build/spectral.txt", "rb") as f:
	for line in f:
		parse = line.strip().split(" ")
		sp_p.append(	float(parse[0]))
		alpha.append(	float(parse[1]))
		nu.append(	-float(parse[2]))
		i+=1
#plt.tricontour(alpha, nu, sp_p, 15, linewidths=0.5, colors='k')
plt.figure(1)
plt.tricontourf(alpha,nu,sp_p, 15,cmap=plt.cm.jet)
plt.colorbar()



plt.figure(2)
xi = np.linspace(0,2,5000)
yi = np.linspace(0,3,5000)
zi = griddata((alpha, nu), sp_p, (xi[None,:], yi[:,None]), method='linear')
plt.contour(xi, yi, zi)
CS = plt.contourf(xi,yi,zi,15,cmap=plt.cm.jet,locator=ticker.LogLocator())
plt.colorbar()


plt.show()
