from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def func(xs, a, b, c, d):
	#return a*(1-xs) + b*(1-xs)**2 + c*(1-xs)**3 + d
	return (a)*((1-xs)**b) + c
	#return a + b*xs + c*np.exp(d*(1-xs));
	#return a*(1-xs) + c*np.exp(d*(1-xs));

x = []; l = []; u = [];
with open("f2nf2p.dat","rb") as f:
	for line in f:
		parse = line.strip().split(" ")
		if( float(parse[0]) > 0.95 ): continue
		x.append( float(parse[0]) )
		l.append( float(parse[1]) )
		u.append( float(parse[2]) )
x = np.asarray(x); l = np.asarray(l); u = np.asarray(u);
m = (l+u)/2.
err = (l-u)/2.
popt, pcov = curve_fit(func, x, m, sigma=1./(err))
print popt

plt.fill_between(x, l, u)
plt.errorbar(x,m,yerr=err,fmt='ro')
x = np.arange(0,1,0.001)
plt.plot(x, func(x, *popt), color='black', linestyle='--',linewidth=3)
plt.show()
