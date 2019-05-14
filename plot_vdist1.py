#Python modules
import h5py
import numpy as np
import matplotlib.pyplot as plt

#from Numeric import array, matrixmultiply, transpose
#from LinearAlgebra import inverse

#from mpi4py import MPI
import struct

# Athena++ modules
import athena_read
#import get_ave
from smooth import *
from sys import byteorder

#np.set_printoptions(threshold=np.nan)

fig = plt.figure()

tstar = 100.0
dim = 2
lz = 1000

cstar = (8.314510E7 * tstar)**0.5
hstar = 3.0857E17
#gstar = 4.19251E-7 * Sigma
gstar = cstar * cstar / hstar

#files=['../kt.out1.00075.athdf','../kt.out1.00150.athdf','../kt.out1.00225.athdf','../kt.out1.00300.athdf']
files=['../kt.out1.00100.athdf','../kt.out1.00200.athdf','../kt.out1.00300.athdf']

Loop = 0
for Loop in range(3):
  filename = files[Loop]

  leveln=2

  quantities=['rho','vel1','vel2','vel3','pgas','Er','Fr1','Fr2','Fr3','Pr11','Pr22','Pr33','Pr12','Pr13','Pr23','Er0','Fr01','Fr02','Fr03']

  f=h5py.File(filename, 'r')

  attributes = f.attrs.items()
  attrs = dict(attributes)
  level = f.attrs['MaxLevel']

  time = f.attrs['Time']
  subsample = False
  if leveln is not None:
     if level > leveln:
      subsample = True
      level = leveln

  #print "Real Athena++ athdf file from level {:d}".format(level)
  data = athena_read.athdf(filename, quantities=quantities, level=level, subsample=subsample)

  f.close()

  dens_tot=data['rho']

  critical = 2.0e-4
  #vmean = -138.19021332199961

  vbkgd = data['vel1'][np.where(dens_tot < critical)]
  vmean = np.mean(vbkgd)

  dens= dens_tot[np.where(dens_tot > critical)]
  vx  = data['vel1'][np.where(dens_tot > critical)]-vmean
  vy  = data['vel2'][np.where(dens_tot > critical)]
  vz  = data['vel3'][np.where(dens_tot > critical)]
  Vc = np.sqrt(vx*vx+vy*vy+vz*vz)
  Mc = dens*Vc

  #print Vc
  #print max(Vc), min(Vc)

  NN = 500

  hist, edges = np.histogram(Vc, range=[-10,800],bins=NN)
  # Histogram
  #print hist
  #print edges

  i=0

  outfile=open('v_dist', 'w')

  total_dens = np.sum(dens)
  #print total_dens

  #p=[]

  for l, r in zip(edges[:-1], edges[1:]):
      den_sum = np.sum(dens[(Vc > l) & (Vc < r)])/total_dens
      #p.append= np.sum(dens[(Vc > l) & (Vc < r)])/total_dens
      #if (i==6):
        #print i,dens[(Vc > l) & (Vc < r)],den_sum
      #print i,den_sum,0.5*(l+r)
      outfile.write(format(i)+'  ')
      outfile.write(format(0.5*(l+r))+'  ')
      outfile.write(format(den_sum)+'  ')
      outfile.write('\n')
      i=i+1

  outfile.close()

  data=np.loadtxt('v_dist')

  if (Loop==0):
   x0=data[:,1]
   y0=data[:,2]
   #sum = 0.0
   #i=0
   #for i in range(NN):
    #sum=sum+data[i,2]
    #print sum
   yhat0=smooth(y0)
  elif (Loop==1):
   x1=data[:,1]
   y1=data[:,2]
   #sum = 0.0
   #for i in range(NN):
    #sum=sum+data[i,2]
    #print sum
   yhat1=smooth(y1)
  elif (Loop==2):
   x2=data[:,1]
   y2=data[:,2]
   #sum = 0.0
   #for i in range(NN):
    #sum=sum+data[i,2]
    #print sum
   yhat2=smooth(y2)
  elif (Loop==3):
   x3=data[:,1]
   y3=data[:,2]
   #sum = 0.0
   #for i in range(NN):
    #sum=sum+data[i,2]
    #print sum
   yhat3=smooth(y3)
 
  Loop =Loop+1

plt.plot(x0,yhat0,'k-',color='black')
plt.plot(x1,yhat1,'k--',color='red')
plt.plot(x2,yhat2,'k:',color='blue')
#plt.plot(x3,yhat3,'k-.')
plt.legend([r"$t=2\,t_*$",r"$t=4\,t_*$",r"$t=6\,t_*$"], handlelength=3.5)
ax = plt.gca()
ax.get_yaxis().set_tick_params(direction='in')
ax.get_xaxis().set_tick_params(direction='in')
plt.xlim([0,700])
plt.ylim([0,0.03])
ax11 = ax.twinx()
ax11.set_ylim(ax.get_ylim())
ax11.get_yaxis().set_tick_params(direction='in',labelright='off')

plt.savefig('dist.pdf', format='pdf')




