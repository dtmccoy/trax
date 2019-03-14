from get_mac import *
from numpy.ma import *
import datetime
from numpy import *
from trax import *
from get_merra2_slp import *
time_in=[2005,1,1]
time_av=get_month_ymd(time_in)
	
SLP0=get_slp(time_in)#[time_in.year,time_in.month,time_in.day])
time_av=get_month_ymd(time_in)
SLPM=zeros((len(time_av),SLP0[0].shape[0],SLP0[0].shape[1]))
for  i in range(len(time_av)):
	SLPM[i]=get_slp(time_av[i])[0]
SLPP=SLP0[0]-nanmean(SLPM,axis=0)
DD=GET_LAT_LON(SLPP,SLP0[0],SLP0[1],SLP0[2])
DD2=GET_LAT_LON(SLP0[0],SLP0[0],SLP0[1],SLP0[2])
LWP=get_mac(time_in,'clwp')
COMPOS=do_composite(LWP[0],LWP[1],LWP[2],[2005,1,1],DD)
CSLP=do_composite(SLP0[0],SLP0[1],SLP0[2],[2005,1,1],DD)
#pcolormesh(SLP0[2],SLP0[1],masked_invalid(SLP0[0]))
pcolormesh(LWP[2],LWP[1],masked_invalid(LWP[0]))
scatter(DD2['lon'],DD2['lat'],marker='x',c='r')
scatter(DD['lon'],DD['lat'],marker='+',c='r')
