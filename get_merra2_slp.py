def get_slp(time):
        import iris
        from netCDF4 import Dataset
        from time2str import *
	import glob
	from numpy import *
	from conv_25slp import *
	wd='/group_workspaces/jasmin2/asci/dtmccoy/MERRA2_METEO/SLP/'
        fn=wd+'MERRA2_*.tavg1_2d_slv_Nx.' + time2str(time) +'.SUB.nc4'
        fn=glob.glob(fn)
        Z=iris.load(fn[0])
        Z=Z[1]
        lat=1.*Z.coord('latitude').points
        lon=1.*Z.coord('longitude').points
        lon[lon>180]=lon[lon>180]-360
        ind=argsort(lon)
        lon=lon[ind]
        Z=squeeze(Z[:,:,ind].data)
        ind2=argsort(lat)
        Z=Z[ind2]
        lat=lat[ind2]
	return Z,lat,lon
