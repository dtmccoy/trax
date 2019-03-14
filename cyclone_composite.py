from scipy.stats import *
from matplotlib.pyplot import *
from numpy import *
from numpy.ma import *
def get_xy_from_latlon(lon,lat):
	import numpy as np
	import copy
	x=np.zeros((lon.shape[0],lat.shape[0]))
	y=copy.copy(x)
	R=6371.0
	#lon[lon>180]=lon[lon>180]-360
	lon=lon/180.0*np.pi
	lat=lat/180.0*np.pi
	lono=copy.copy(x)
	lato=copy.copy(x)
	for i in range(0,lon.shape[0]):
		x[i,:]=lon[i]
		lono[i,:]=180*lon[i]/np.pi
	for i in range(0,lat.shape[0]):
		y[:,i]=lat[i]
		lato[:,i]=180*lat[i]/np.pi
	xd=copy.copy(x)
	yd=copy.copy(y)
	xd=np.cos(y)*R*x
	yd=y*R
	return (xd,yd,lato,lono)
	
def locate_cyclone_centers(x,y,lon,lat,slp,crit,slpp):
	import copy
	import numpy as np
	#timedim=1
	#if len(slp.shape)>2:
	timedim=slp.shape[2]
	latcent=[]
	loncent=[]
	#if timedim>1:
	for i in range(0,timedim):
		latcentt,loncentt=locate_centers_snapshot_curve(x,y,lat,lon,slp[:,:,i],crit,slpp[:,:,i])
		latcent.append(latcentt)
		loncent.append(loncentt)
	"""
	else:
		latcentt,loncentt=locate_centers_snapshot_curve(x,y,lat,lon,slp,crit)
                latcent.append(latcentt)
                loncent.append(loncentt)
	"""
	return latcent,loncent

def locate_centers_snapshot_curve(x,y,lat,lon,slp,crit,slpp):
	import numpy as np
	import copy
	dx=np.gradient(x,axis=0)
	dy=np.gradient(y,axis=1)
	dslp=np.gradient(slpp)
	dslp[0]=dslp[0]/dx
	dslp[1]=dslp[1]/dy
	ddslp=copy.copy(dslp)
	ddslp[0]=np.gradient(dslp[0],axis=0)
	ddslp[1]=np.gradient(dslp[1],axis=1)
	ddslp[0]=ddslp[0]/dx
	ddslp[1]=ddslp[1]/dy
	dxdyslp=np.gradient(dslp[0],axis=1)/dy
	dxdyslp=masked_invalid(dxdyslp)
	ddslp=masked_invalid(ddslp)
	dslp=masked_invalid(dslp)
	slp=masked_invalid(slp)
	'''
	from matplotlib.pyplot import *
	figure()
	pcolormesh(masked_invalid(ddslp[0]+ddslp[1]))
	colorbar()
	figure()
	pcolormesh(masked_invalid(dxdyslp))
	colorbar()
	figure()
	pcolormesh(masked_invalid(masked_invalid(slp)))
	colorbar()
        figure()
        pcolormesh(masked_invalid(masked_invalid(slpp)))
        colorbar()
	'''
	ind_candidate=(np.abs(lat)>30)&(np.abs(lat)<80) & (slp<crit[0]) \
	& ((ddslp[0]+ddslp[1])>crit[1]) & (dxdyslp<crit[2])
	#ind_candidate=(np.abs(lat)>30)&(np.abs(lat)<80) & (slp<crit[0])&((ddslp[0]+ddslp[1])<crit[1]) & (dxdyslp>crit[2])
	lonind=lon[ind_candidate]
	latind=lat[ind_candidate]	
	slpind=slp[ind_candidate]
	latcent, loncent=get_apart(latind,lonind,slpind,crit)
	
	return latcent, loncent	
def get_apart(latind,lonind,slpind,crit):
        good=ones(lonind.shape,dtype=bool)
	for i in range(0,lonind.size):
                dist=get_dd(latind,lonind,latind[i],lonind[i])#RE*(arccos(sin(lonindr)*sin(lonindr[i])+cos(lonindr)*cos(lonindr[i])*cos(latindr-latindr[i])))
                closemin=dist<crit[3]
                MINSLP=np.amin(slpind[closemin])
                good_ss=good[closemin]
                good_ss[slpind[closemin]>MINSLP]=False
                good[closemin]=good_ss
	latcent=latind[good]
        loncent=lonind[good]
	return latcent, loncent
def get_dd(latx,lonx,lat0x,lon0x):
	lat=latx/180.*pi; lon=lonx/180.*pi; lat0=lat0x/180.*pi; lon0=lon0x/180.*pi
	RE=6378.
	dist=RE*(arccos(sin(lon)*sin(lon0)+cos(lon)*cos(lon0)*cos(lat-lat0)))	
	return dist
	
def composite_plan_centers(lon,lat,latcent,loncent,z,time,stat='mean'):
	import numpy as np
	lens=[18,19]
	#lens=[30,31]
	ncy=len(loncent)
	sz=2000
	plangrid=(np.linspace(-sz,sz,lens[0]),np.linspace(-sz,sz,lens[1]))
	comp=np.zeros([lens[0]-1,lens[1]-1,ncy])
	num=np.zeros([lens[0]-1,lens[1]-1,ncy])
	comp[:]=np.nan
	tick=0
	statsc=np.zeros((ncy,9))*np.NaN
	for jj in range(0,len(loncent)):
		TT=time
		comp[:,:,jj],num[:,:,jj],statstemp=comp_centers2(lon,lat,latcent[jj],loncent[jj], z,comp[:,:,jj],num[:,:,jj],plangrid,TT,stat=stat)
		statsc[tick,:]=statstemp
		tick=tick+1
	return comp,num,plangrid,statsc
def comp_centers2(lon,lat,latcent,loncent,z,comp,num,plangrid,time,stat='mean'):	
	if loncent>180:
		loncent=loncent-360
	lon[lon>180]=lon[lon>180]-360
	R=6371.0
        #loncentn=lon[np.argmin(np.abs(lon[:,0]-loncent)),0]
        #latcentn=lat[0,np.argmin(np.abs(lat[0,:]-latcent))]
	dlon=lon-loncent
        dlat=lat-latcent
        dlon[dlon>180]=dlon[dlon>180]-360 # make sure that center includes the globe
        dlon[dlon<-180]=dlon[dlon<-180]+360
	dy=(dlat/180)*pi*R
	dx=cos(lat/180.*pi)*R*(dlon/180.*pi)
	ind=isnan(z)==False
	z=z[ind]
	dx=dx[ind]
	dy=dy[ind]
	statso=binned_statistic_2d(dx,dy,z,bins=[plangrid[0],plangrid[1]],statistic=stat);
	statsn=binned_statistic_2d(dx,dy,z,bins=(plangrid[0],plangrid[1]),statistic='count');
	ZZ=z[sqrt(dx**2+dy**2)<2000.]
	NUM=len(where(isnan(ZZ)==False))
	ind=isnan(ZZ)==False
	statscyc=[np.mean(ZZ[ind]),NaN,NaN,NUM,latcent,loncent,time[0],time[1],time[2]]
	return statso[0],statsn[0],statscyc
				
