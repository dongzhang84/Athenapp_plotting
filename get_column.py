def get_column(filename,outfile,tstar,dim,lz):
    # Python modules
    import h5py
    import numpy as np

    #from mpi4py import MPI
    import struct

    # Athena++ modules
    import athena_read
    #import get_ave
    from sys import byteorder

    tstar = 1000.0
    dim = 2
    lz = 3000
  
    cstar = (8.314510E7 * tstar)**0.5
    hstar = 3.0857E17
    #gstar = 4.19251E-7 * Sigma
    gstar = cstar * cstar / hstar

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

    # Determine new grid size
    nx1 = attrs['RootGridSize'][0] * 2**level
    nx2 = attrs['RootGridSize'][1] * 2**level
    nx3 = attrs['RootGridSize'][2] * 2**level

    f.close()

    dens_tot=data['rho']
    #tgas_tot=data['pgas']/dens
    #vx_tot=data['vel1'][np.where(dens0 <2.0e-4)]
    #vy_tot=data['vel2']
    #vz_tot=data['vel3']

    i=0
    column =0

    while (i<=999):
     if (dens_tot[0,199,i]>1.99e-4):
        #print column,i
        column = column + dens_tot[0,199,i]
     i=i+1   

    critical = 2.0e-4
    #vmean = -138.19021332199961

    dens= dens_tot[np.where(dens_tot > critical)]


    table=[time,column]

    print(table)
    outfile.write(format(time)+'  ')
    outfile.write(format(column)+'  ')
    outfile.write('\n')

    return


