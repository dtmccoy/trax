def get_mac(time,nm):
        from lon180 import *
	import iris
        wd='/group_workspaces/jasmin2/asci/dtmccoy/MACLWP_daily/'
        monstr=str(time[1])
        if time[1]<10:
                monstr='0'+monstr
        fn=wd+nm+'1deg_maclwpv1.'+str(time[0])+monstr+'.nc4'#201606.nc4
        Z=iris.load(fn)
        INDER=0
        #print Z
        if nm=='tlwp':
                INDER=3
        Z=Z[INDER]#[time[2]]
        lat=1.*Z.coord('latitude').points
        lon=1.*Z.coord('longitude').points
        #fn=nm+'.npz'
        #fh=load(wd+fn)
        #lat=fh['lat']
        #lon=fh['lon']
        #timef=fh['time']
        #ind=(time[0]==timef[:,0])&(time[1]==timef[:,1])&(time[2]==timef[:,2])
        #data=squeeze(fh['data'][ind])
        lon,data=lon180(lon,Z[time[2]-1].data)
        data[data<0]=NaN
        return data,lat,lon
