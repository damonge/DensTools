import numpy as np
import sys as sys
import matplotlib.pyplot as plt
from matplotlib import cm

def read_grid(prefix,nfiles) :
    f=open(prefix+".0000","rb")
    num_grids,ngrid=np.fromfile(f,dtype=np.int32,count=2)
    f.close()

    print "Will read %d fields"%num_grids+" with %d^3 nodes"%ngrid
    grid_out=np.zeros([ngrid,ngrid,ngrid,num_grids])
    for ifil in np.arange(nfiles) :
        f=open(prefix+".%04d"%ifil,"rb")
        nug,ng=np.fromfile(f,dtype=np.int32,count=2)
        if (nug!=num_grids) or (ng!=ngrid) :
            print "shit"
            sys.exit(1)
        nx_here=np.fromfile(f,dtype=np.int32,count=1)[0]
        print "File #%d"%(ifil+1)+", %d slices found"%nx_here
        for ix in np.arange(nx_here) :
            ix_this=np.fromfile(f,dtype=np.int32,count=1)[0]
            grid_out[ix_this,:,:,:]=np.fromfile(f,dtype=np.float32,count=ng*ng*nug).reshape([ng,ng,nug])
        f.close()

    if num_grids==1 :
        grid_out=grid_out[:,:,:,0]

    return grid_out

dens=read_grid("output_dens",2)
dens_sm=read_grid("output_dens_sm",2)
vel=read_grid("output_vel",2)
lvel=read_grid("output_linvel",2)
tid=read_grid("output_tidal",2)
ng=len(dens)
i_slice=ng/2

fig,ax=plt.subplots()
cax=ax.imshow(dens[i_slice,:,:],interpolation='none',origin='lower')
cbar=fig.colorbar(cax)

fig,ax=plt.subplots()
cax=ax.imshow(dens_sm[i_slice,:,:],interpolation='none',origin='lower')
cbar=fig.colorbar(cax)

fig,ax=plt.subplots()
cax=ax.imshow(tid[i_slice,:,:,0]+tid[i_slice,:,:,1]+tid[i_slice,:,:,2],interpolation='none',origin='lower')
cbar=fig.colorbar(cax)

fig,ax=plt.subplots()
cax=ax.imshow(np.sqrt(vel[i_slice,:,:,0]**2+vel[i_slice,:,:,1]**2+vel[i_slice,:,:,2]**2),interpolation='none',origin='lower')
cbar=fig.colorbar(cax)

fig,ax=plt.subplots()
cax=ax.imshow(np.sqrt(lvel[i_slice,:,:,0]**2+lvel[i_slice,:,:,1]**2+lvel[i_slice,:,:,2]**2),interpolation='none',origin='lower')
cbar=fig.colorbar(cax)
plt.show()
