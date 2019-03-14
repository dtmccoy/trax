def conv_25slp_fr(SLP,latSLP,lonSLP,latd,lond,griddataf=False):
        from numpy import *
        from scipy.interpolate import interp2d
	if len(where(isnan(SLP))[0])>0:
		print 'detected NAN'
		griddataf=True # otherwise interp2d breaks
	if griddataf==False:
                SLPi=interp2d(lonSLP,latSLP,SLP)#,array(lond),array(latd),)
                SLPdt=SLPi(array(lond),array(latd))
	        SLPi=interp2d(array(lond),array(latd),SLPdt) # add this section to deal with the fact that FW07 uses 2.5deg data interpolated to 1deg
                SLPd=SLPi(lonSLP,latSLP)
        else:
                #latd=arange(-90,90,2.5); lond=arange(-180,180,2.5)
                #bb=(latd,lond)
                #from scipy.stats import binned_statistic_dd
                #datao2=binned_statistic_dd((lat,lon),datao,bins=bb)
                from scipy.interpolate import griddata
                [xx,yy]=meshgrid(lonSLP,latSLP)
                z=ravel(SLP)
                xx=ravel(xx)
                yy=ravel(yy)
                ind=isnan(z)==False
                SLPd= griddata((xx[ind], yy[ind]), z[ind], (lond[None,:], latd[:,None]), method='linear')
		SLPi=interp2d(lond,latd,SLPd)
		SLPd=SLPi(lonSLP,latSLP)
	latSLP=latd; lonSLP=lond
        return SLPd,latSLP,lonSLP
