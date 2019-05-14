# Python modules
import h5py
import numpy as np

#from mpi4py import MPI
import struct

# Athena++ modules
#import athena_read
#import get_ave
from get_func import *
from sys import byteorder

tstar = 100.0
dim = 2
nx = 1000
ny = 400
nz= 1
lz = nx

lab='../kt.out1'
#outfile = 'cloud.column'
outfile = 'cloud.ave'

#outfile=open('cloud.column', 'w')
outfile=open('cloud.ave', 'w')

i=0

while (i<=400):
  if (i>=0 and i<10):
    filename=lab+'.0000'+str(i)+'.athdf'
  elif (i>=10 and i<100):
    filename=lab+'.000'+str(i)+'.athdf'
  elif (i>=100 and i<1000):
    filename=lab+'.00'+str(i)+'.athdf'
  print filename
  #print get_column2(filename,outfile,tstar,dim,nx,ny,nz,lz)
  print get_ave(filename,outfile,tstar,dim,lz)
  i=i+1

outfile.close()

