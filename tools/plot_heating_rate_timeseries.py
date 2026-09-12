import numpy as np
import matplotlib.pyplot as plt
expid='REF010'
za,heatrate=np.loadtxt('../data_out/radiative_heating_rate_exp_'+expid+'.dat',unpack=True)
dt=60.0 # model time step
nlev=80
nstep=int(48*3600/dt) + 1
x1 = np.arange(0,nstep)
y1 = np.arange(0,nlev)
X1, Y1 = np.meshgrid(x1,y1)

plt.xlabel('Model time steps')
plt.ylabel('Level number')
plt.ylim(nlev-1,50,20)
plt.xlim(0,nstep-1)
plt.xticks(np.arange(0,nstep,400))
plt.yticks(np.arange(0,nlev-1,5))
xvar=np.reshape(heatrate,(nstep,nlev))
var1=xvar.transpose()
plt.title('Longwave heating rate (K/day) - '+expid)
figure = plt.gcf()
figure.set_size_inches(7, 4)
plt.contourf(X1,Y1,var1,20,cmap='jet_r')
plt.colorbar()
plt.contour(X1,Y1,var1,20,colors='black',linewidths=0.5)
plt.savefig('../plots/heating_rate_'+expid+'.png',dpi=600)
plt.show()


