import numpy as np
import matplotlib.pyplot as plt
expid='REF015'
dt=900.0
nlevs=14
nstep=int(48*3600/dt) + 1
zs,ts,ws,ls,dw=np.loadtxt('../data_out/soil_profiles_exp_'+expid+'.dat',unpack=True)
xvar=np.nstep=int(48*3600/dt) + 1
xvar=np.reshape(ts,(nstep,nlevs))
var1=xvar.transpose()
x1 = np.arange(0,nstep)
y1 = np.arange(0,nlevs)
X1, Y1 = np.meshgrid(x1,y1)

plt.xlabel('Time steps')
plt.ylabel('Level number')
plt.ylim(nlevs-1,0,-1)
plt.xlim(0,nstep-1)
plt.xticks(np.arange(0,nstep,24))
plt.yticks(np.arange(0,nlevs-1,1))
xvar=np.reshape(ts,(nstep,nlevs))
var1=xvar.transpose()
plt.title('Soil temperature (K) - '+expid)
figure = plt.gcf()
figure.set_size_inches(7, 4)
plt.contourf(X1,Y1,var1,20,cmap='coolwarm')
plt.colorbar()
plt.contour(X1,Y1,var1,20,colors='black',linewidths=0.5)
plt.savefig('../plots/soil_temperature_'+expid+'.png',dpi=600)
plt.show()


plt.xlabel('Time steps')
plt.ylabel('Level number')
plt.ylim(nlevs-1,0,-1)
plt.xlim(0,nstep-1)
plt.xticks(np.arange(0,nstep,24))
plt.yticks(np.arange(0,nlevs-1,1))
xvar=np.reshape(ws,(nstep,nlevs))
var1=xvar.transpose()
plt.title('Soil moisture content (m$^3$/m$^3$) - '+expid)
figure = plt.gcf()
figure.set_size_inches(7, 4)
plt.contourf(X1,Y1,var1,20,cmap='Spectral')
plt.colorbar()
plt.contour(X1,Y1,var1,20,colors='black',linewidths=0.5)
plt.savefig('../plots/soil_moisture_'+expid+'.png',dpi=600)
plt.show()

plt.xlabel('Time steps')
plt.ylabel('Level number')
plt.ylim(nlevs-1,0,-1)
plt.xlim(0,nstep-1)
plt.xticks(np.arange(0,nstep,24))
plt.yticks(np.arange(0,nlevs-1,1))
xvar=np.reshape(ls,(nstep,nlevs))
var1=xvar.transpose()
plt.title('Thermal conductivity (m$^2$s$^{-1}$) - '+expid)
figure = plt.gcf()
figure.set_size_inches(7, 4)
plt.contourf(X1,Y1,var1,20,cmap='Spectral')
plt.colorbar()
plt.contour(X1,Y1,var1,20,colors='black',linewidths=0.5)
plt.savefig('../plots/soil_thermal_conductivity_'+expid+'.png',dpi=600)
plt.show()


