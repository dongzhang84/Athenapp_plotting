def get_ave(filename,outfile,tstar,dim,lz):
    # Python modules
    import h5py
    import numpy as np

    #from mpi4py import MPI
    import struct

    # Athena++ modules
    import athena_read
    from sys import byteorder

    tstar = tstar
    dim = dim
    lz = lz

    #filename = '../kt.out1.00100.athdf'
    #outfile = 'hotwind.ave'

    #print filename
  
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

    critical = 2.0e-4
    #vmean = -138.19021332199961

    dens= dens_tot[np.where(dens_tot > critical)]
    tgas= data['pgas'][np.where(dens_tot > critical)]/dens
    vx  = data['vel1'][np.where(dens_tot > critical)]
    vy  = data['vel2'][np.where(dens_tot > critical)]
    vz  = data['vel3'][np.where(dens_tot > critical)]

    vbkgd= data['vel1'][np.where(dens_tot <= critical)]
    vmean= np.mean(vbkgd)

    #print vx

    mass=np.mean(dens)
    #totmass=7868.3543189564225
    totmass=np.sum(dens)
    cloud = np.sum(dens_tot[np.where(dens_tot >= 1.0/3)])


    if (dim == 2):
       vx0=np.mean(vx*dens)/mass
       vy0=np.mean(vy*dens)/mass
       sigma0=np.sqrt(np.mean(dens*((vx-vx0)**2.0))/mass)
       sigma1=np.sqrt(np.mean(dens*((vy-vy0)**2.0))/mass)
       sigma =np.sqrt(sigma0**2.0+sigma1**2.0)
       fcloud=cloud/totmass
    else:
       vx0=np.mean(vx*dens)/mass
       vy0=np.mean(vy*dens)/mass
       vz0=np.mean(vz*dens)/mass
       sigma0=np.sqrt(np.mean(dens*((vx-vx0)**2.0))/mass)
       sigma1=np.sqrt(np.mean(dens*((vy-vy0)**2.0))/mass)
       sigma2=np.sqrt(np.mean(dens*((vz-vz0)**2.0))/mass)
       sigma =np.sqrt(sigma0**2.0+sigma1**2.0+sigma2**2.0)
       fcloud=cloud/totmass

    table=[time, mass, vx0, vy0, sigma0, sigma1, sigma, fcloud,vmean]

    print(table)
    #print(outfile)
    outfile.write(format(time)+'  ')
    outfile.write(format(mass)+'  ')
    outfile.write(format(vx0)+'  ')
    outfile.write(format(vy0)+'  ')
    outfile.write(format(sigma0)+'  ')
    outfile.write(format(sigma1)+'  ')
    outfile.write(format(sigma)+'  ')
    outfile.write(format(fcloud)+'  ')
    outfile.write(format(vmean)+'  ')
    outfile.write('\n')

    return


##=============================================================================##


def get_column(filename,outfile,tstar,dim,nx,ny,nz,lz):
    # Python modules
    import h5py
    import numpy as np

    #from mpi4py import MPI
    import struct

    # Athena++ modules
    import athena_read
    #import get_ave
    from sys import byteorder

    tstar = tstar
    dim = dim
    nx = nx
    ny = ny
    nz = nz
    lz = lz
  
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
    ny_half = ny/2-1

    nz_half = 0
    if (dim == 3):
      nz_half = nz/2-1

    while (i<=nx-1):
     if (dens_tot[nz_half,ny_half-1,i]>1.99e-4):
        #print column,i
        column = column + dens_tot[nz_half,ny_half-1,i]
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




##=============================================================================##


def get_column2(filename,outfile,tstar,dim,nx,ny,nz,lz):
    # Python modules
    import h5py
    import numpy as np

    #from mpi4py import MPI
    import struct

    # Athena++ modules
    import athena_read
    #import get_ave
    from sys import byteorder

    tstar = tstar
    dim = dim
    nx = nx
    ny = ny
    nz = nz
    lz = lz
  
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
    nx_half = nx/2-1
    
    nz_half = 0
    if (dim == 3):
      nz_half = nz/2-1

    while (i<=ny-1):
     if (dens_tot[nz_half,i,nx_half]>1.99e-4):
        #print column,i
        column = column + dens_tot[nz_half,i,nx_half]
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
