from lon180 import *
from numpy import *
import cyclone_composite
## read in data from some variable and composite it on SLP. 
## SLP and Z are in lon,lat
## return composites,num in each bin, the grid, statistics on mean values, and where the lat and lon of each composite are centered.
def do_composite(Z,lat,lon,time,CENTERS,stat='mean'):
	#x,y,lat_g,lon_g=cyclone_composite.get_xy_from_latlon(lon,lat)
        lat_g,lon_g=meshgrid(lat,lon,indexing='ij')
	compos,num,plangrid,stats=cyclone_composite.composite_plan_centers(lon_g,lat_g,CENTERS['lat'],CENTERS['lon'],Z,array(time),stat=stat)
	return compos,num,plangrid,stats,CENTERS['lat'],CENTERS['lon']
def GET_LAT_LON(SLP,SLP0,latslp,lonslp):
	from conv_25slp_fr import *
        latd=arange(-90,90,2.5); lond=arange(-180,180,2.5)
	SLP=lon180(lonslp,SLP)[1]
	lonslp,SLP0=lon180(lonslp,SLP0)
	
	SLP=conv_25slp_fr(SLP,latslp,lonslp,latd,lond)[0]	
	if len(where(isnan(SLP0))[0])>0: #deal with NANS in SLP
                from scipy.interpolate import griddata
                [xx,yy]=meshgrid(lonslp,latslp)
                z=ravel(SLP0)
                xx=ravel(xx)
                yy=ravel(yy)
                ind=isnan(z)==False
                SLP0= griddata((xx[ind], yy[ind]), z[ind], (lonslp[None,:], latslp[:,None]), method='linear')	
	from get_land_mask import *
	LSMASK=get_land_mask(latslp,lonslp)
	SLP=SLP*LSMASK
	SLP0=SLP0*LSMASK
	if SLP.shape[0]!=len(lonslp):
                SLP=transpose(SLP)
		SLP0=transpose(SLP0)
	slp=SLP
	slp=expand_dims(slp,axis=2);
        SLP0=expand_dims(SLP0,axis=2);
        x,y,lat_g,lon_g=cyclone_composite.get_xy_from_latlon(lonslp,latslp)
	crit=[101500,6e-3,3e-3,2000]
	latcent,loncent=cyclone_composite.locate_cyclone_centers(x,y,lon_g,lat_g,slp,crit,SLP0)
	return {'lat':latcent[0],'lon':loncent[0]}
def get_month_ymd(time_in,calendar='standard',units='days since 0001-01-01 00:00:00'):
	import datetime
	from netcdftime import utime
	cdftime = utime(units,calendar=calendar)
	time_in=datetime.datetime(int(time_in[0]),int(time_in[1]),int(time_in[2]))
	dt=arange(-15,16)
	time_av=zeros((len(dt),3))
        NUM=cdftime.date2num(time_in)
	for i in range(len(dt)):
		to=cdftime.num2date(NUM+(dt[i]))#,calendar=calendar,units=units)
		time_av[i]=[to.year,to.month,to.day]
	return time_av
