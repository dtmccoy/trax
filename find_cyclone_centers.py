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
def locate_cyclone_centers(x,y,lon,lat,slp,crit):
        import copy
        import numpy as np
        #timedim=1
        #if len(slp.shape)>2:
        timedim=slp.shape[2]
        latcent=[]
        loncent=[]
        #if timedim>1:
        for i in range(0,timedim):
                latcentt,loncentt=locate_centers_snapshot_curve(x,y,lat,lon,slp[:,:,i],crit)
                latcent.append(latcentt)
                loncent.append(loncentt)
        return latcent,loncent
def locate_centers_snapshot_curve(x,y,lat,lon,slp,crit):
        slp[slp<700]=NaN
        import numpy as np
        import copy
        dx=np.gradient(x,axis=0)
        dy=np.gradient(y,axis=1)
        dslp=np.gradient(slp)
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
        ind_candidate=(np.abs(lat)>30)&(np.abs(lat)<80) & (slp<crit[0]) \
        & ((ddslp[0]+ddslp[1])>crit[1]) & (dxdyslp<crit[2])
        xind=x[ind_candidate]
        yind=y[ind_candidate]
        lonind=lon[ind_candidate]
        latind=lat[ind_candidate]
        slpind=slp[ind_candidate]
        ## go through and make sure far enough apart
        RE=6378.
        lonindr=lonind/180.*pi
        latindr=latind/180.*pi
        for i in range(0,lonind.size):
                dist=RE*(arccos(sin(lonindr)*sin(lonindr[i])+cos(lonindr)*cos(lonindr[i])*cos(latindr-latindr[i])))
                subind=dist<crit[3]
                if len(slpind[subind])>0:
                        if slpind[i]==np.amin(slpind[subind]):
                                slpind[i]=-1
        latcent=latind[slpind==-1]
        loncent=lonind[slpind==-1]
        return latcent, loncent
def seg_centers(loncent,latcent,latrange):
        import copy
        loncent1=copy.copy(loncent)
        latcent1=copy.copy(latcent)
        loncent2=copy.copy(loncent)
        latcent2=copy.copy(latcent)
        import numpy as np
        for i in range(0,len(loncent)):
                ind=np.where((latcent[i]>np.min(latrange)) & (latcent[i]<np.max(latrange)))
                latcent1[i]=latcent[i][ind]
                loncent1[i]=loncent[i][ind]
                latcent2[i]=np.delete(latcent[i],ind)
                loncent2[i]=np.delete(loncent[i],ind)
        return latcent2,loncent2,latcent1,loncent1
